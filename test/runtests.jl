using AbstractStorage, AzSessions, AzStorage, CloudSeis, Dates, Distributed, FolderStorage, HTTP, JSON, Random, Test, UUIDs

# trigger extensions
using Blosc
if !Sys.iswindows()
    using CvxCompress
end
using ZfpCompression

global session
if haskey(ENV, "CLIENT_ID") && haskey(ENV, "CLIENT_SECRET") && haskey(ENV, "CLIENT_ID")
    AzSessions.write_manifest(;client_id=ENV["CLIENT_ID"], client_secret=ENV["CLIENT_SECRET"], tenant=ENV["TENANT_ID"])
    session = AzSession(;protocal=AzClientCredentials, resource="https://storage.azure.com/")
else
    session = AzSession(;scope="openid+offline_access+https://storage.azure.com/user_impersonation")
end

const storageaccount1 = ENV["STORAGE_ACCOUNT1"]
const storageaccount2 = ENV["STORAGE_ACCOUNT2"]

abstract type Cloud end
struct Azure <: Cloud end
struct Azure2 <: Cloud end
struct POSIX <: Cloud end

mkcontainer(::Type{Azure}, foldername) = AzContainer(foldername; session=session, storageaccount=storageaccount1)
mkcontainer(::Type{Azure2}, foldername) = [AzContainer(foldername; session=session, storageaccount=a) for a in (storageaccount1, storageaccount2)]
mkcontainer(::Type{POSIX}, foldername) = Folder(foldername)

function csopen_robust(containers, mode; kwargs...)
    local io
    timeout = 180
    tic = time()
    while true
        try
            if (mode == "r" || mode == "r+") && isa(containers, AbstractArray)
                io = csopen(containers[1], mode; kwargs...)
            else
                io = csopen(containers, mode; kwargs...)
            end
            break
        catch e
            @warn "caught exception in csopen, sleeping for 60 seconds"
            io = IOBuffer()
            showerror(io, e)
            write(io, "\n\terror type: $(typeof(e))")
            if VERSION >= v"1.7"
                for (exc, bt) in current_exceptions(current_task())
                    showerror(io, exc, bt)
                    println(io)
                end
            end
            s = String(take!(io))
            close(io)
            @warn s
            sleep(60)
            if time() - tic > timeout
                throw(e)
            end
        end
    end
    io
end

const clouds = (Azure, Azure2, POSIX)
const compressors = Sys.iswindows() ? ("none","blosc","leftjustify","zfp") : ("none","blosc","leftjustify","zfp","cvx")

@testset "CloudSeis, cloud=$cloud, compresser=$compressor" for cloud in clouds, compressor in compressors
    if compressor == "cvx"
        compressor_options = (bx=16,by=16,bz=16,scale=0.0f0) # accurate as possible for cvx unit tests
    end
    compressor_options = () # lossless compression for zfp unit-tests

    rtol = compressor == "cvx" ? 1e-2 : 1e-7

    @testset "CloudSeis, compression selection" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor, compressor_options=compressor_options)
        if compressor == "none"
            @test isa(io.cache.compressor, CloudSeis.NotACompressor)
            @test Dict(io.cache.compressor)["method"] == "none"

            # default compression is "leftjustify"
            r = uuid4()
            io_default = csopen_robust(mkcontainer(cloud, "test-$r-default-cs"), "w", axis_lengths=[10,11,12])
            @test isa(io_default.cache.compressor, CloudSeis.LeftJustifyCompressor)
            rm(io_default)
        end
        if compressor == "blosc"
            @test isa(io.cache.compressor, CloudSeis.BloscCompressor)
            @test Dict(io.cache.compressor)["method"] == "blosc"
        end
        if compressor == "zfp"
            @test isa(io.cache.compressor, CloudSeis.ZfpCompressor)
            @test Dict(io.cache.compressor)["method"] == "zfp"
        end
        if compressor == "cvx"
            @test isa(io.cache.compressor, CloudSeis.CloudSeisCvxCompressor)
            @test Dict(io.cache.compressor)["method"] == "cvx"
        end
        if compressor == "leftjustify"
            @test isa(io.cache.compressor, CloudSeis.LeftJustifyCompressor)
            @test Dict(io.cache.compressor)["method"] == "leftjustify"
        end
        rm(io)
    end

    @testset "allocframe" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor, compressor_options=compressor_options)
        t,h = allocframe(io)

        @test size(t) == (10,11)
        @test eltype(t) == Float32

        @test size(h) == (io.hdrlength,11)
        @test eltype(h) == UInt8

        close(io)
        rm(io)
    end

    @testset "writeframe/readframe with headers and variable fold" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor, compressor_options=compressor_options)

        x = rand(Float32,10,11,12)
        t,h = allocframe(io)

        n = size(h,2)
        J = Vector{Vector{Int}}(undef, size(io,3))
        for iframe = 1:size(io,3)
            t .= x[:,:,iframe]
            if iframe == 2
                J[iframe] = [1:11;]
            elseif iframe == 5
                J[iframe] = Int[]
            elseif iframe == 9
                J[iframe] = [6]
            else
                J[iframe] = randperm(n)[1:div(n,2)]
            end
            for itrace = 1:size(io,2)
                set!(prop(io,io.axis_propdefs[2]), h, itrace, itrace)
                set!(prop(io,io.axis_propdefs[3]), h, itrace, iframe)
                set!(prop(io,stockprop[:TRC_TYPE]), h, itrace, itrace ∈ J[iframe] ? tracetype[:live] : tracetype[:dead])
                if itrace ∉ J[iframe]
                    t[:,itrace] .= 0
                end
            end
            writeframe(io,t,h)
        end
        close(io)

        nbytes = filesize(io.extents[1].container, io.extents[1].name)
        if compressor == "leftjustify"
            @test nbytes == size(io,3)*8 + mapreduce(iframe->length(J[iframe])*(size(io,1)*sizeof(io.traceformat) + headerlength(io)), +, 1:size(io,3))
        elseif compressor == "none"
            @test nbytes == size(io,3)*8 + size(io,3)*size(io,2)*(4*size(io,1) + headerlength(io))
        elseif compressor == "blosc"
            @test nbytes < size(io,3)*8 + size(io,3)*size(io,2)*(4*size(io,1) + headerlength(io))
        elseif compressor == "zfp"
            @test nbytes < size(io,3)*8 + size(io,3)*size(io,2)*(4*size(io,1) + headerlength(io))
        end

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        for iframe = 1:size(io,3)
            t,h = readframe(io, iframe)
            for i = 1:size(h,2)
                if i ∈ J[iframe]
                    @test t[:,i] ≈ x[:,i,iframe] rtol=rtol
                    @test get(prop(io,io.axis_propdefs[3]), h, i) == iframe
                    @test get(prop(io,io.axis_propdefs[2]), h, i) == i
                end
                @test get(prop(io,stockprop[:TRC_TYPE]), h, i) == (i ∈ J[iframe] ? tracetype[:live] : tracetype[:dead])
            end
        end
        rm(io)
    end

    @testset "writeframe/readframe with index" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], dataproperties=[DataProperty("P",1)], compressor=compressor, compressor_options=compressor_options)

        x = rand(Float32,10,11)

        t,h = allocframe(io)
        t .= x

        writeframe(io, t, 1)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        t,h = readframe(io, 1)
        close(io)
        @test t ≈ x rtol=rtol
        rm(io)
    end

    @testset "write, lstarts=$lstarts, lincs=$lincs" for (lstarts,lincs) = (((1,1,1),(1,1,1)), ((2,3,4),(5,6,7)))
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lstarts=lstarts, axis_lincs=lincs, compressor=compressor, compressor_options=compressor_options)

        x = rand(Float32,10,11,12)
        for ifrm = 1:12
            x[:,:,ifrm] .= ifrm
        end
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        for ifrm = 1:12
            _x = readframetrcs(io, lstarts[3] + lincs[3]*(ifrm - 1))
            @test _x ≈ x[:,:,ifrm] rtol=rtol
        end
        close(io)
        rm(io)
    end

    @testset "write, multiple extents" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], frames_per_extent=1, compressor=compressor, compressor_options=compressor_options)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        for ifrm = 1:12
            _x = readframetrcs(io, ifrm)
            @test _x ≈ x[:,:,ifrm] rtol=rtol
        end
        close(io)
        rm(io)
    end

    @testset "readtrcs, lstarts=$lstarts, lincs=$lincs" for (lstarts,lincs) = (((1,1,1),(1,1,1)), ((2,3,4),(5,6,7)))
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lstarts=lstarts, axis_lincs=lincs, compressor=compressor, compressor_options=compressor_options)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        _x = readtrcs(io, :, :, :)
        close(io)
        rm(io)

        @test x ≈ _x rtol=rtol
    end

    @testset "readhdrs, lstarts=$lstarts, lincs=$lincs" for (lstarts,lincs) = (((1,1,1),(1,1,1)), ((2,3,4),(5,6,7)))
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lincs=lincs, axis_lstarts=lstarts, compressor=compressor, compressor_options=compressor_options)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        h = readhdrs(io, :, :, :)
        close(io)
        rm(io)

        for ifrm = 1:12, itrc = 1:11
            @test get(prop(io, "TRACE"), view(h, :, :, ifrm), itrc) ≈ lstarts[2] + lincs[2]*(itrc - 1)
            @test get(prop(io, "FRAME"), view(h, :, :, ifrm), itrc) ≈ lstarts[3] + lincs[3]*(ifrm - 1)
        end
    end

    @testset "readhdrs!, lstarts=$lstarts, lincs=$lincs" for (lstarts,lincs) = (((1,1,1),(1,1,1)), ((2,3,4),(5,6,7)))
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lincs=lincs, axis_lstarts=lstarts, compressor=compressor, compressor_options=compressor_options)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        h = readhdrs(io, :, :, :)
        close(io)
        rm(io)

        for ifrm = 1:12, itrc = 1:11
            @test get(prop(io, "TRACE"), view(h, :, :, ifrm), itrc) ≈ lstarts[2] + lincs[2]*(itrc - 1)
            @test get(prop(io, "FRAME"), view(h, :, :, ifrm), itrc) ≈ lstarts[3] + lincs[3]*(ifrm - 1)
        end
    end

    @testset "partial read for foldmap, lstarts=$lstarts, lincs=$lincs" for (lstrts,lncs) = (((1,1,1),(1,1,1)), ((2,3,4),(5,6,7)))
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lstarts=lstrts, axis_lincs=lncs, compressor=compressor, compressor_options=compressor_options)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        f = [fold(io,i) for i in lstarts(io,3) .+ lincs(io,3)*(0:size(io,3)-1)]
        close(io)
        rm(io)

        @test f == [11 for i=1:12]
    end

    @testset "full read for foldmap, lstarts=$lstarts, lincs=$lincs" for (lstrts,lncs) = (((1,1,1),(1,1,1)), ((2,3,4),(5,6,7)))
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lstarts=lstrts, axis_lincs=lncs, compressor=compressor, compressor_options=compressor_options, frames_per_extent=2)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        f = foldmap(io; all=true)
        @test f == [11 for i=1:12]
    end

    @testset "similarto" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w",
            axis_pincs=[0.1,0.2,0.3],
            axis_lengths=[10,11,12],
            dataproperties=[DataProperty("P",1)],
            compressor=compressor,
            compressor_options=compressor_options)
        close(io)

        _io = csopen_robust(mkcontainer(cloud, "test-$r-sim-cs"), "w", similarto=mkcontainer(cloud, "test-$r-cs"))
        @test size(_io) == (10,11,12)
        @test pincs(_io)[1] ≈ 0.1
        @test pincs(_io)[2] ≈ 0.2
        @test pincs(_io)[3] ≈ 0.3

        rm(io)
        rm(_io)
    end

    @testset "similarto, changing extents" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w",
            axis_pincs=[0.1,0.2,0.3],
            axis_lengths=[10,11,12],
            dataproperties=[DataProperty("P",1)],
            compressor=compressor,
            compressor_options=compressor_options)
        @test length(io.extents) == 1
        close(io)

        _io = csopen_robust(mkcontainer(cloud, "test-$r-sim-cs"), "w", similarto=mkcontainer(cloud, "test-$r-cs"), frames_per_extent=1)
        @test size(_io) == (10,11,12)
        @test pincs(_io)[1] ≈ 0.1
        @test pincs(_io)[2] ≈ 0.2
        @test pincs(_io)[3] ≈ 0.3

        @test length(_io.extents) == 12 
    end

    @testset "propdefs" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor, compressor_options=compressor_options)
        @test propdefs(io,1).label == "SAMPLE"
        @test propdefs(io,2).label == "TRACE"
        @test propdefs(io,3).label == "FRAME"
        @test propdefs(io)[1].label == "SAMPLE"
        @test propdefs(io)[2].label == "TRACE"
        @test propdefs(io)[3].label == "FRAME"
        close(io)
        rm(io)
    end

    @testset "props" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor, compressor_options=compressor_options)
        @test props(io,1).def.label == "SAMPLE"
        @test props(io,2).def.label == "TRACE"
        @test props(io,3).def.label == "FRAME"
        @test props(io)[1].def.label == "SAMPLE"
        @test props(io)[2].def.label == "TRACE"
        @test props(io)[3].def.label == "FRAME"
        close(io)
        rm(io)
    end

    @testset "vector props" begin
        r = uuid4()
        pdef = TracePropertyDef("X","XX",Vector{Float64},2)
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], tracepropertydefs=[pdef], compressor=compressor, compressor_options=compressor_options)
        t,h = allocframe(io)
        for i = 1:11
            set!(prop(io,"TRACE"), h, i, i)
            set!(prop(io,"FRAME"), h, i, 1)
            set!(prop(io,"TRC_TYPE"), h, i, tracetype[:live])
            set!(prop(io,"X"), h, i, i*[1.0,2.0])
        end
        writeframe(io,t,h)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        t,h = readframe(io,1)

        for i = 1:11
            @test get(prop(io,"X"), h, i) ≈ i*[1.0,2.0]
        end
        close(io)
        rm(io)
    end

    @testset "string props" begin
        r = uuid4()
        pdef = TracePropertyDef("X","XX",Vector{UInt8},32)
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], tracepropertydefs=[pdef], compressor=compressor, compressor_options=compressor_options)
        t,h = allocframe(io)
        for i = 1:11
            set!(prop(io,"TRACE"), h, i, i)
            set!(prop(io,"FRAME"), h, i, 1)
            set!(prop(io,"TRC_TYPE"), h, i, tracetype[:live])
            set!(prop(io,"X"), h, i, "HELLO")
        end
        writeframe(io,t,h)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        t,h = readframe(io,1)

        for i = 1:11
            @test unsafe_string(pointer(get(prop(io,"X"), h, i))) == "HELLO"
        end
        close(io)
        rm(io)
    end

    @testset "pincs" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_pincs=[1.0,2.0,3.0], compressor=compressor, compressor_options=compressor_options)
        @test pincs(io,1) ≈ 1.0
        @test pincs(io,2) ≈ 2.0
        @test pincs(io,3) ≈ 3.0
        @test pincs(io)[1] ≈ 1.0
        @test pincs(io)[2] ≈ 2.0
        @test pincs(io)[3] ≈ 3.0
        close(io)
        rm(io)
    end

    @testset "pstarts" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_pstarts=[1.0,2.0,3.0], compressor=compressor, compressor_options=compressor_options)
        @test pstarts(io,1) ≈ 1.0
        @test pstarts(io,2) ≈ 2.0
        @test pstarts(io,3) ≈ 3.0
        @test pstarts(io)[1] ≈ 1.0
        @test pstarts(io)[2] ≈ 2.0
        @test pstarts(io)[3] ≈ 3.0
        close(io)
        rm(io)
    end

    @testset "lincs" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lincs=[1,2,3], compressor=compressor, compressor_options=compressor_options)
        @test lincs(io,1) == 1
        @test lincs(io,2) == 2
        @test lincs(io,3) == 3
        @test lincs(io)[1] == 1
        @test lincs(io)[2] == 2
        @test lincs(io)[3] == 3
        close(io)
        rm(io)
    end

    @testset "lstarts" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lstarts=[1,2,3], compressor=compressor, compressor_options=compressor_options)
        @test lstarts(io,1) == 1
        @test lstarts(io,2) == 2
        @test lstarts(io,3) == 3
        @test lstarts(io)[1] == 1
        @test lstarts(io)[2] == 2
        @test lstarts(io)[3] == 3
        close(io)
        rm(io)
    end

    @testset "units" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_units=["X","Y","Z"], compressor=compressor, compressor_options=compressor_options)
        @test units(io,1) == "X"
        @test units(io,2) == "Y"
        @test units(io,3) == "Z"
        @test units(io)[1] == "X"
        @test units(io)[2] == "Y"
        @test units(io)[3] == "Z"
        close(io)
        rm(io)
    end

    @testset "domains" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_domains=["X","Y","Z"], compressor=compressor, compressor_options=compressor_options)
        @test domains(io,1) == "X"
        @test domains(io,2) == "Y"
        @test domains(io,3) == "Z"
        @test domains(io)[1] == "X"
        @test domains(io)[2] == "Y"
        @test domains(io)[3] == "Z"
        close(io)
        rm(io)
    end

    @testset "azimuth information" begin 
        x_direction = "east"
        y_direction = "north"
        z_direction = "elevation"
        tti_angle_units = "degrees"
        tti_azimuth_positive_direction = "clockwise"
        tti_azimuth_origin_axis = "y"
        tti_symmetry_axis_z_direction = "elevation"

        
        g = Geometry(
            ox=1.0,oy=2.0,oz=3.0,
            ux=4.0,uy=5.0,uz=6.0,
            vx=7.0,vy=8.0,vz=9.0,
            wx=10.0,wy=11.0,wz=12.0,
            u1=1,un=2,
            v1=1,vn=4,
            w1=1,wn=6,
            umin = 0.0, umax = 10.0,
            vmin = 1.0, vmax = 11.0,
            wmin = 2.0, wmax = 12.0,
            x_direction=x_direction,
            y_direction=y_direction,
            z_direction=z_direction,
            tti_angle_units = tti_angle_units,
            tti_azimuth_positive_direction = tti_azimuth_positive_direction,
            tti_azimuth_origin_axis = tti_azimuth_origin_axis,
            tti_symmetry_axis_z_direction=tti_symmetry_axis_z_direction)


        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w",axis_propdefs=[stockprop[:IU],stockprop[:IV],stockprop[:IW],stockprop[:VOLUME]],axis_lengths=[g.un,g.vn,g.wn], geometry=g, compressor=compressor, compressor_options=compressor_options)

        _g = geometry(io)
        @test _g.ox ≈ g.ox
        @test _g.oy ≈ g.oy
        @test _g.oz ≈ g.oz
        @test _g.ux ≈ g.ux
        @test _g.uy ≈ g.uy
        @test _g.uz ≈ g.uz
        @test _g.vx ≈ g.vx
        @test _g.vy ≈ g.vy
        @test _g.vz ≈ g.vz
        @test _g.wx ≈ g.wx
        @test _g.wy ≈ g.wy
        @test _g.wz ≈ g.wz
        @test _g.u1 === g.u1
        @test _g.un === g.un
        @test _g.v1 === g.v1
        @test _g.vn === g.vn
        @test _g.w1 === g.w1
        @test _g.wn === g.wn
        @test _g.umin ≈ 0.0
        @test _g.umax ≈ 10.0
        @test _g.vmin ≈ 1.0
        @test _g.vmax ≈ 11.0
        @test _g.wmin ≈ 2.0
        @test _g.wmax ≈ 12.0
        @test _g.x_direction == g.x_direction
        @test _g.y_direction == g.y_direction
        @test _g.z_direction == g.z_direction
        @test _g.tti_angle_units == tti_angle_units
        @test _g.modelstride_u == 1
        @test _g.modelstride_v == g.un
        @test _g.modelstride_w == g.vn*g.un
        @test _g.tti_azimuth_positive_direction == g.tti_azimuth_positive_direction
        @test _g.tti_azimuth_origin_axis == g.tti_azimuth_origin_axis
        @test _g.tti_symmetry_axis_z_direction == g.tti_symmetry_axis_z_direction
        close(io)
        rm(io)
    end

    @testset "geometry backwards compatability" begin
        g = Geometry(
            ox=1.0,oy=2.0,oz=3.0,
            ux=4.0,uy=5.0,uz=6.0,
            vx=7.0,vy=8.0,vz=9.0,
            wx=10.0,wy=11.0,wz=12.0,
            u1=1,un=2,
            v1=1,vn=4,
            w1=1,wn=6)

            r = uuid4()
            io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w",axis_propdefs=[stockprop[:IU],stockprop[:IV],stockprop[:IW],stockprop[:VOLUME]],axis_lengths=[g.un,g.vn,g.wn], geometry=g, compressor=compressor, compressor_options=compressor_options)

            _g = geometry(io)

            @test _g.ox ≈ g.ox
            @test _g.oy ≈ g.oy
            @test _g.oz ≈ g.oz
            @test _g.ux ≈ g.ux
            @test _g.uy ≈ g.uy
            @test _g.uz ≈ g.uz
            @test _g.vx ≈ g.vx
            @test _g.vy ≈ g.vy
            @test _g.vz ≈ g.vz
            @test _g.wx ≈ g.wx
            @test _g.wy ≈ g.wy
            @test _g.wz ≈ g.wz
            @test _g.u1 === g.u1
            @test _g.un === g.un
            @test _g.v1 === g.v1
            @test _g.vn === g.vn
            @test _g.w1 === g.w1
            @test _g.wn === g.wn
            @test _g.umin ≈ 0
            @test _g.umax ≈ 1
            @test _g.vmin ≈ 0
            @test _g.vmax ≈ 1
            @test _g.wmin ≈ 0
            @test _g.wmax ≈ 1
            close(io)
            rm(io)

    end

    @testset "in" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor, compressor_options=compressor_options)

        @test in(stockprop[:TRACE], io) == true
        @test in(stockprop[:CDP],io) == false

        close(io)
        rm(io)
    end

    @testset "dataproperty" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], dataproperties=[DataProperty("X",1)], compressor=compressor, compressor_options=compressor_options)

        @test dataproperty(io, "X") == 1
        @test dataproperty(io, "X", 2) == 1
        @test dataproperty(io, "Y", 2) == 2

        close(io)
        rm(io)
    end

    @testset "hasdataproperty" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], dataproperties=[DataProperty("X",1)], compressor=compressor, compressor_options=compressor_options)

        @test hasdataproperty(io, "X") == true
        @test hasdataproperty(io, "Y") == false

        close(io)
        rm(io)
    end

    @testset "backward compatibility, dataproperty (1.2->1.1)" begin
        r = uuid4()
        c = mkcontainer(cloud, "test-$r-cs")
        _c = isa(c, Array) ? c[1] : c

        io = csopen_robust(c, "w", axis_lengths=[10,11,12], dataproperties=[DataProperty("X",1)], compressor=compressor, compressor_options=compressor_options)
        d = JSON.parse(read(_c, "description.json", String))
        @test isa(d["dataproperties"], Dict)
        d["dataproperties"] = [Dict("label"=>"X", "value"=>1)]
        write(_c, "description.json", json(d, 1))
        close(io)

        io = csopen_robust(c, "r")
        d = JSON.parse(read(_c, "description.json", String))
        @test isa(d["dataproperties"], Vector)
        @test dataproperty(io, "X") == 1

        close(io)
        rm(io)
    end

    @testset "copy!" begin
        r1 = uuid4()
        io1 = csopen_robust(mkcontainer(cloud, "test-$r1-cs"), "w", axis_lengths=[10,11,12], compressor=compressor, compressor_options=compressor_options)
        r2 = uuid4()
        io2 = csopen_robust(mkcontainer(cloud, "test-$r2-cs"), "w", axis_lengths=[10,11,12], compressor=compressor, compressor_options=compressor_options)

        h1 = allocframehdrs(io1)
        h2 = allocframehdrs(io2)

        for i = 1:11
            set!(prop(io1,stockprop[:TRACE]), h1, i, i)
            set!(prop(io1,stockprop[:FRAME]), h1, i, 2)
            set!(prop(io1,stockprop[:TRC_TYPE]), h1, i, tracetype[:live])
        end

        copy!(io2, h2, io1, h1)

        for i = 1:11
            @test get(prop(io2,stockprop[:TRACE]), h2, i) == i
            @test get(prop(io2,stockprop[:FRAME]), h2, i) == 2
            @test get(prop(io2,stockprop[:TRC_TYPE]), h2, i) == tracetype[:live]
        end

        close(io1)
        rm(io1)
        close(io2)
        rm(io2)
    end

    @testset "reduce" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,50], frames_per_extent=2, compressor=compressor, compressor_options=compressor_options)

        @test length(io.extents) == 25

        T = rand(Float32,10,11,50)
        H = rand(UInt8,headerlength(io),11,50)
        for i = 1:50
            writeframe(io, T[:,:,i], i)
            readframehdrs!(io, view(H,:,:,i), i)
        end
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        reduce(io)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        for i = 1:50
            t,h = readframe(io, i)
            @test t ≈ T[:,:,i] rtol=rtol
            @test h == H[:,:,i]
        end

        @test length(io.extents) == 1

        rm(io)
    end

    @testset "reduce, parallel" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,50], frames_per_extent=2, compressor=compressor, compressor_options=compressor_options)

        @test length(io.extents) == 25

        T = rand(Float32,10,11,50)
        H = rand(UInt8,headerlength(io),11,50)
        for i = 1:50
            writeframe(io, T[:,:,i], i)
            readframehdrs!(io, view(H,:,:,i), i)
        end
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        addprocs(5)
        @everywhere using AzStorage, Blosc, CloudSeis, FolderStorage, ZfpCompression

        if !Sys.iswindows()
            @everywhere using CvxCompress
        end

        reduce(io; frames_per_extent=10)
        rmprocs(workers())
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        for i = 1:50
            t,h = readframe(io, i)
            @test t ≈ T[:,:,i] rtol=rtol
            @test h == H[:,:,i]
        end

        @test length(io.extents) == 5

        rm(io)
    end

    @testset "cp" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,50], frames_per_extent=10, compressor=compressor, compressor_options=compressor_options)
        x = rand(Float32,10,11,50)
        write(io, x, :, :, :)
        close(io)
        r = lowercase(randstring('a':'z',6))

        c = mkcontainer(cloud, "test-$r-cs")
        cp(io, c)
        rm(io)

        iocp = csopen(c)
        @test readtrcs(iocp, :, :, :) ≈ x rtol=rtol
        rm(iocp)
    end

    @testset "cp, parallel" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,50], frames_per_extent=10, compressor=compressor, compressor_options=compressor_options)
        x = rand(Float32,10,11,50)
        write(io, x, :, :, :)
        close(io)
        
        r = uuid4()
        c = mkcontainer(cloud, "test-$r-cs")

        addprocs(2)
        @everywhere using AzStorage, Blosc, CloudSeis, FolderStorage, ZfpCompression

        if !Sys.iswindows()
            @everywhere using CvxCompress
        end

        cp(io, c)
        rmprocs(workers())

        rm(io)

        iocp = csopen(c)
        @test readtrcs(iocp, :, :, :) ≈ x rtol=rtol
        rm(iocp)
    end

    @testset "mv" begin
        r = uuid4()
        c_in = mkcontainer(cloud, "test-$r-cs")
        io = csopen_robust(c_in, "w", axis_lengths=[10,11,50], frames_per_extent=10, compressor=compressor, compressor_options=compressor_options)
        x = rand(Float32,10,11,50)
        write(io, x, :, :, :)
        close(io)

        r = uuid4()
        c = mkcontainer(cloud, "test-$r-cs")

        mv(io, c)

        iomv = csopen(c)
        @test readtrcs(iomv, :, :, :) ≈ x rtol=rtol
        if isa(c_in, AbstractArray)
            @test isdir(c_in[1]) == false
        else
            @test isdir(c_in) == false
        end
        rm(iomv)
    end

    @testset "robust cscreate" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], force=true, compressor=compressor, compressor_options=compressor_options)
        writeframe(io, rand(Float32,10,11), 1)
        close(io)
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], force=true, compressor=compressor, compressor_options=compressor_options)
        writeframe(io, rand(Float32,10,11), 1)
        close(io)
    end

    @testset "backward compatibility (1.1->1.0)" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], force=true, compressor="none", compressor_options=compressor_options)
        description = JSON.parse(read(io.containers[1], "description.json", String))
        delete!(description, "compressor")
        write(io.containers[1], "description.json", json(description))
        io = csopen(mkcontainer(cloud, "test-$r-cs"), axis_lengths=[10,11,12], force=true, compressor="none", compressor_options=compressor_options)
        @test isa(io.cache, CloudSeis.Cache{CloudSeis.NotACompressor})
    end

    @testset "similar with new axis lengths" begin
        r = uuid4()
        container = mkcontainer(cloud, "test-$r-cs")
        io = csopen_robust(container, "w", axis_lengths=[10,11,2])

        trcs1 = rand(Float32,10,11)
        trcs2 = rand(Float32,10,11)

        writeframe(io, trcs1, 1)
        writeframe(io, trcs2, 2)
        close(io)

        container_similar = mkcontainer(cloud, "test-$r-similar-cs")
        io_similar = csopen_robust(container_similar, "w"; similarto=container, axis_lengths=[10,11,1])
        writeframe(io_similar, trcs1, 1)
        close(io_similar)

        io = csopen_robust(container, "r")
        io_similar = csopen_robust(container_similar, "r")
        _trcs1 = readframetrcs(io_similar, 1)
        @test _trcs1 ≈ trcs1 rtol=rtol

        rm(io)
        rm(io_similar)
    end

    @testset "similarto with new number of dimensions" begin
        r = uuid4()

        container = mkcontainer(cloud, "test-$r-cs")

        io = csopen_robust(container, "w",
            axis_pincs=[0.1,0.2,0.3,0.4],
            axis_lengths=[10,11,12,13],
            dataproperties=[DataProperty("P",1)],
            compressor=compressor)

        io2 = csopen_robust(mkcontainer(cloud, "test-$r-sim-cs"), "w", similarto=container, axis_lengths=[10,11,12])
        description = JSON.parse(read(io2.containers[1], "description.json", String))
        @test length(description["fileproperties"]["axis_lengths"]) == 3
        @test length(description["fileproperties"]["axis_propdefs"]) == 3
        @test length(description["fileproperties"]["axis_units"]) == 3
        @test length(description["fileproperties"]["axis_domains"]) == 3
        @test length(description["fileproperties"]["axis_pstarts"]) == 3
        @test length(description["fileproperties"]["axis_pincs"]) == 3
        @test length(description["fileproperties"]["axis_lstarts"]) == 3
        @test length(description["fileproperties"]["axis_lincs"]) == 3

        @test lincs(io2) == (1,1,1)
        @test lstarts(io2) == (1,1,1)
        @test pincs(io2) == (0.1,0.2,0.3)
        @test pstarts(io2) == (0.0,0.0,0.0)
        @test domains(io2) == (stockdomain[:UNKNOWN],stockdomain[:UNKNOWN],stockdomain[:UNKNOWN])
        @test units(io2) == (stockunit[:UNKNOWN],stockunit[:UNKNOWN],stockunit[:UNKNOWN])
        @test propdefs(io2) == (SAMPLE=stockprop[:SAMPLE], TRACE=stockprop[:TRACE], FRAME=stockprop[:FRAME])

        @test_throws ErrorException csopen(mkcontainer(cloud, "test-$r-sim-2-cs"), "w", similarto=container, axis_lengths=[10,11,12,13,14])

        rm(io)
        rm(io2)
    end

    @testset "csopen for 6 dimensional data-set" begin
        r = uuid4()
        container = mkcontainer(cloud, "test-$r-cs")
        io = csopen_robust(container, "w", axis_lengths=[10,11,2,2,3,1])
        close(io)
        io = csopen_robust(container, "r")
        @test size(io) == (10,11,2,2,3,1)
        labels = [propdefs(io)[i].label for i=1:6]
        @test labels == ["SAMPLE", "TRACE", "FRAME", "VOLUME", "HYPRCUBE", "DIM6"]
    end

    @testset "non unitary logical start and increments" begin
        r = uuid4()
        container = mkcontainer(cloud, "test-$r-cs")

        io = csopen_robust(container, "w", axis_lengths=[12,11,4], axis_lstarts=[1,2,5], axis_lincs=[1,1,3])

        # test map from cartesian to linear index
        @test CloudSeis.linearframeidx(io, 5) == 1
        @test CloudSeis.linearframeidx(io, 8) == 2
        @test CloudSeis.linearframeidx(io, 11) == 3
        @test CloudSeis.linearframeidx(io, 14) == 4

        # test map from unitary to non-unitary cartesian index
        @test CloudSeis.logicalframeidx(io, 1) == CartesianIndex((5,))
        @test CloudSeis.logicalframeidx(io, 2) == CartesianIndex((8,))
        @test CloudSeis.logicalframeidx(io, 3) == CartesianIndex((11,))
        @test CloudSeis.logicalframeidx(io, 4) == CartesianIndex((14,))

        @test CloudSeis.lineartraceidx(io, 2) == 1
        @test CloudSeis.logicaltraceidx(io, 1) == 2

        rm(io)
    end

    @testset "non unitary logical start and increments, read/write frame" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lstarts=[0,2,6], axis_lincs=[1,2,3], compressor=compressor, compressor_options=compressor_options)

        x = rand(Float32,10,11,12)
        t,h = allocframe(io)

        n = size(h,2)
        J = Vector{Vector{Int}}(undef, size(io,3))
        for iframe = 1:size(io,3)
            t .= x[:,:,iframe]
            if iframe == 2
                J[iframe] = [1:11;]
            elseif iframe == 5
                J[iframe] = Int[]
            elseif iframe == 9
                J[iframe] = [6]
            else
                J[iframe] = randperm(n)[1:div(n,2)]
            end
            for itrace = 1:size(io,2)
                set!(prop(io,io.axis_propdefs[2]), h, itrace, 2 + 2*(itrace - 1))
                set!(prop(io,io.axis_propdefs[3]), h, itrace, 6 + 3*(iframe - 1))
                set!(prop(io,stockprop[:TRC_TYPE]), h, itrace, itrace ∈ J[iframe] ? tracetype[:live] : tracetype[:dead])
                if itrace ∉ J[iframe]
                    t[:,itrace] .= 0
                end
            end
            writeframe(io,t,h)
        end
        close(io)

        nbytes = filesize(io.extents[1].container, io.extents[1].name)
        if compressor == "leftjustify"
            @test nbytes == size(io,3)*8 + mapreduce(iframe->length(J[iframe])*(size(io,1)*sizeof(io.traceformat) + headerlength(io)), +, 1:size(io,3))
        elseif compressor == "none"
            @test nbytes == size(io,3)*8 + size(io,3)*size(io,2)*(4*size(io,1) + headerlength(io))
        elseif compressor == "blosc"
            @test nbytes < size(io,3)*8 + size(io,3)*size(io,2)*(4*size(io,1) + headerlength(io))
        elseif compressor == "zfp"
            @test nbytes < size(io,3)*8 + size(io,3)*size(io,2)*(4*size(io,1) + headerlength(io))
        elseif compressor == "cvx"
            @test nbytes < size(io,3)*8 + size(io,3)*size(io,2)*(4*size(io,1) + headerlength(io))
        end

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        for iframe = 1:size(io,3)
            t,h = readframe(io, 6 + 3*(iframe-1))
            for i = 1:size(h,2)
                if i ∈ J[iframe]
                    @test t[:,i] ≈ x[:,i,iframe] rtol=rtol
                    @test get(prop(io,io.axis_propdefs[3]), h, i) == 6 + 3*(iframe-1)
                    @test get(prop(io,io.axis_propdefs[2]), h, i) == 2 + 2*(i - 1)
                end
                @test get(prop(io,stockprop[:TRC_TYPE]), h, i) == (i ∈ J[iframe] ? tracetype[:live] : tracetype[:dead])
            end
        end
        rm(io)
    end

    @testset "non unitary logical start and increments, writeframe/readframe with index" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lstarts=[0,2,6], axis_lincs=[1,2,3], dataproperties=[DataProperty("P",1)], compressor=compressor, compressor_options=compressor_options)

        x = rand(Float32,10,11)

        t,h = allocframe(io)
        t .= x

        writeframe(io, t, 9)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        t,h = readframe(io, 9)
        close(io)
        @test t ≈ x rtol=rtol
        rm(io)
    end

    @testset "non unitary logical start and increments, exception for in-between out-of-bounds index" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lstarts=[1,2,3], axis_lincs=[4,5,6], dataproperties=[DataProperty("P",1)], compressor=compressor, compressor_options=compressor_options)

        @test_throws ErrorException CloudSeis.linearframeidx(io, 4)
        @test @inbounds CloudSeis.linearframeidx(io, 2) == 1

        @test_throws ErrorException CloudSeis.lineartraceidx(io, 3)
        @test @inbounds CloudSeis.linearframeidx(io, 3) == 1
    end

    @testset "backwards compatability, logical increments and starts" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], force=true)
        description = JSON.parse(read(io.containers[1], "description.json", String))
        delete!(description, "axis_lincs")
        delete!(description, "axis_lstarts")
        write(io.containers[1], "description.json", json(description))
        io = csopen(mkcontainer(cloud, "test-$r-cs"), axis_lengths=[10,11,12], force=true)
        @test lstarts(io) == (1,1,1)
        @test lincs(io) == (1,1,1)
    end

    @testset "frame iterator, LogicalIndices" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], force=true)
        idxs = CartesianIndex[]
        for idx in LogicalIndices(io)
            push!(idxs, idx)
        end
        _idxs = [CartesianIndex((i,)) for i=1:12]
        @test idxs == _idxs
        rm(io)

        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_lstarts=[1,1,3], axis_lincs=[1,1,4], force=true)
        idxs = CartesianIndex[]
        for idx in LogicalIndices(io)
            push!(idxs, idx)
        end
        _idxs = [CartesianIndex((3 + (i-1)*4,)) for i=1:12]
        @test idxs == _idxs
        rm(io)

        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,4,6], force=true)
        idxs = CartesianIndex[]
        for idx in LogicalIndices(io)
            push!(idxs, idx)
        end
        _idxs = [CartesianIndex((i,j)) for i=1:4, j=1:6][:]
        @test idxs == _idxs
        rm(io)

        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,4,6], axis_lstarts=[1,1,2,3], axis_lincs=[1,1,4,5], force=true)
        idxs = CartesianIndex[]
        for idx in LogicalIndices(io)
            push!(idxs, idx)
        end
        _idxs = [CartesianIndex((2+(i-1)*4,3+(j-1)*5)) for i=1:4, j=1:6][:]
        @test idxs == _idxs
        rm(io)
    end

    @testset "partial read for foldmap with empty extents" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], force=true)
        @test fold(io,1) == 0
    end

    @testset "left justify headers" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,12,1], force=true)

        h = allocframehdrs(io)
        for itrace = 1:12
            set!(prop(io, stockprop[:TRACE]), h, itrace, itrace)
            set!(prop(io, stockprop[:FRAME]), h, itrace, 1)
            if rem(itrace,2) == 0
                set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:live])
            else
                set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:dead])
            end
        end

        CloudSeis.leftjustify!(io, h)
        for itrace = 1:6
            @test get(prop(io, stockprop[:TRC_TYPE]), h, itrace) == tracetype[:live]
            @test get(prop(io, stockprop[:TRACE]), h, itrace) == 2*itrace
        end
        for itrace = 7:12
            @test get(prop(io, stockprop[:TRC_TYPE]), h, itrace) == tracetype[:dead]
        end
    end
end

@testset "CloudSeis optimized read, cloud=$cloud" for cloud in clouds
    r = uuid4()
    io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,12,3], force=true)

    t,h = allocframe(io)

    for iframe = 1:3
        for itrace = 1:12
            set!(prop(io, stockprop[:TRACE]), h, itrace, itrace)
            set!(prop(io, stockprop[:FRAME]), h, itrace, iframe)
            if rem(itrace,2) == 0
                set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:live])
                t[:,itrace] .= itrace*iframe
            else
                set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:dead])
                t[:,itrace] .= 0
            end
        end

        writeframe(io, t, h)
    end

    close(io)
    io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "r")

    for iframe = 1:3
        _t,_h = readframe(io, iframe; regularize=false)

        for itrace = 1:6
            @test _t[:,itrace] ≈ 2*itrace*iframe*ones(Float32,10)
            @test get(prop(io, stockprop[:TRACE]), _h, itrace) == 2*itrace
            @test get(prop(io, stockprop[:FRAME]), _h, itrace) == iframe
            @test get(prop(io, stockprop[:TRC_TYPE]), _h, itrace) == tracetype[:live]
        end

        for itrace = 7:12
            @test get(prop(io, stockprop[:TRC_TYPE]), _h, itrace) == tracetype[:dead]
        end
    end

    rm(io)
end

@testset "CloudSeis optimized read empty frame, cloud=$cloud" for cloud in clouds
    r = uuid4()
    io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,12,2], force=true)

    t,h = allocframe(io)

    for itrace = 1:12
        set!(prop(io, stockprop[:TRACE]), h, itrace, itrace)
        set!(prop(io, stockprop[:FRAME]), h, itrace, 1)
        if rem(itrace,2) == 0
            set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:live])
            t[:,itrace] .= itrace
        else
            set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:dead])
            t[:,itrace] .= 0
        end
    end

    writeframe(io, t, h)

    close(io)
    io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "r")

    _t,_h = readframe(io, 2; regularize=false)

    for itrace = 1:12
        @test get(prop(io, stockprop[:TRC_TYPE]), _h, itrace) == tracetype[:dead]
    end

    rm(io)
end

@testset "CloudSeis optimized read, mixing regularize on/off, cloud=$cloud" for cloud in clouds
    r = uuid4()
    io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,12,2], force=true)

    t,h = allocframe(io)

    for iframe = 1:2
        for itrace = 1:12
            set!(prop(io, stockprop[:TRACE]), h, itrace, itrace)
            set!(prop(io, stockprop[:FRAME]), h, itrace, iframe)
            if rem(itrace,2) == 0
                set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:live])
                t[:,itrace] .= itrace*iframe
            else
                set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:dead])
                t[:,itrace] .= 0
            end
        end

        writeframe(io, t, h)
    end

    close(io)
    io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "r")

    _t,_h = readframe(io, 1; regularize=false)
    _t,_h = readframe(io, 2)

    for itrace = 1:12
        if rem(itrace,2) == 0
            @test _t[:,itrace] ≈ itrace*2*ones(Float32,10)
            @test get(prop(io, stockprop[:TRACE]), _h, itrace) == itrace
            @test get(prop(io, stockprop[:FRAME]), _h, itrace) == 2
            @test get(prop(io, stockprop[:TRC_TYPE]), _h, itrace) == tracetype[:live]
        else
            @test get(prop(io, stockprop[:TRC_TYPE]), _h, itrace) == tracetype[:dead]
        end
    end

    rm(io)

    @testset "left justify" begin
        r = uuid4()
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,12,1], force=true)

        t,h = allocframe(io)
        j = 0
        for itrace = 1:12
            set!(prop(io, stockprop[:TRACE]), h, itrace, itrace)
            set!(prop(io, stockprop[:FRAME]), h, itrace, 1)
            if rem(itrace,2) == 0
                set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:live])
                j += 1
                t[:,itrace] .= j
            else
                set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:dead])
            end
        end
        writeframe(io, t, h)

        CloudSeis.leftjustify!(io, t, h)
        for i = 1:6
            @test t[:,i] ≈ i*ones(Float32, 10)
        end

        CloudSeis.regularize!(io, t, h)

        CloudSeis.leftjustify!(io, t, h, 6)
        for i = 1:6
            @test t[:,i] ≈ i*ones(Float32, 10)
        end

        close(io)
        rm(io)
    end
end

@testset "CloudSeis optimized read only supported for 'r'" for cloud in clouds
    r = uuid4()
    container = mkcontainer(cloud, "test-$r-cs")
    io = csopen_robust(container, "w", axis_lengths=[10,12,1], force=true)

    t,h = allocframe(io)
    j = 0
    for itrace = 1:12
        set!(prop(io, stockprop[:TRACE]), h, itrace, itrace)
        set!(prop(io, stockprop[:FRAME]), h, itrace, 1)
        if rem(itrace,2) == 0
            set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:live])
            j += 1
            t[:,itrace] .= j
        else
            set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:dead])
        end
    end
    writeframe(io, t, h)
    close(io)

    io = csopen(container, "r+")
    @test_throws ErrorException readframe(io, 1; regularize=false)

    rm(io)
end

@testset "Behave well when writing an empty frame" for cloud in clouds
    r = uuid4()
    container = mkcontainer(cloud, "test-$r-cs")
    io = csopen_robust(container, "w", axis_lengths=[10,12,1], force=true)

    t,h = allocframe(io)
    for itrace = 1:size(io,2)
        set!(prop(io, stockprop[:TRC_TYPE]), h, itrace, tracetype[:dead])
    end
    @test writeframe(io, t, h) == 0
    close(io)

    rm(io)
end

@testset "mutate the description of a CloudSeis dataset" for cloud in clouds
    container = mkcontainer(cloud, "test-$(uuid4())-cs")
    io = csopen_robust(container, "w", axis_lengths=[10,12,20,2], frames_per_extent=10, dataproperties=[DataProperty("v",1), DataProperty("ρ",2)])
    n1,n2 = 10,12
    for i4 = 1:2, i3 = 1:20
        writeframe(io, i4*i3*ones(n1,n2), i3, i4)
    end
    close(io)

    io = csopen_robust(container, "r+")

    @test_throws ErrorException description!(io, axis_lengths=[10,12,21,4])

    description!(io, axis_lengths=[10,12,20,4], dataproperties=[DataProperty("v",1), DataProperty("ρ",2), DataProperty("vs",3), DataProperty("epsilon",4)])
    close(io)

    io = csopen_robust(container, "r")
    @test size(io) == (10,12,20,4)
    @test dataproperty(io, "v") == 1
    @test dataproperty(io, "ρ") == 2
    @test dataproperty(io, "vs") == 3
    @test dataproperty(io, "epsilon") == 4
    for i4 = 1:2, i3 = 1:20
        t = readframetrcs(io, i3, i4)
        @test t ≈ i4*i3*ones(n1,n2)
    end
    for i4 = 3:4, i3=1:20
        @test fold(io, i3, i4) == 0
    end

    io = csopen_robust(container, "w", axis_lengths=[10,12,20,2], frames_per_extent=10, dataproperties=[DataProperty("v",1), DataProperty("ρ",2)])
    n1,n2 = 10,12
    for i4 = 1:2, i3 = 1:20
        writeframe(io, i4*i3*ones(n1,n2), i3, i4)
    end
    close(io)

    io = csopen_robust(container, "r+")
    description!(io, axis_lengths=[10,12,20,4], dataproperties_add=[DataProperty("vs",3), DataProperty("epsilon",4)])

    io = csopen_robust(container, "r")
    @test size(io) == (10,12,20,4)
    @test dataproperty(io, "v") == 1
    @test dataproperty(io, "ρ") == 2
    @test dataproperty(io, "vs") == 3
    @test dataproperty(io, "epsilon") == 4

    rm(io)
end

@testset "history construction and access" for cloud in clouds
    x = history!(
        process = "myprocess1",
        process_parameters = Dict("one"=>1,"two"=>2),
        flow_parameters = Dict("three"=>3),
        histories=[])

    x = history!(x, process="myprocess2", process_parameters = Dict("four"=>4,"five"=>5))

    @test !haskey(x, "inputs")
    @test x["flow"]["processes"][1]["process"] == "myprocess1"
    @test x["flow"]["processes"][1]["parameters"] == Dict("one"=>1,"two"=>2)
    @test x["flow"]["processes"][2]["process"] == "myprocess2"
    @test x["flow"]["processes"][2]["parameters"] == Dict("four"=>4,"five"=>5)
    @test x["flow"]["parameters"] == Dict("three"=>3)

    container = mkcontainer(cloud, "test-$(uuid4())-cs")
    io = csopen_robust(container, "w", axis_lengths=[10,12,20,2], frames_per_extent=10, history=x)
    close(io)

    y = history(io)
    @test x == y

    z = history!(
        process = "myprocess3",
        process_parameters = Dict("six"=>6, "seven"=>7),
        histories = [container]
    )

    @test z["inputs"][1]["history"] == y

    if isa(container, Vector)
        @test z["inputs"][1]["container"] == merge(minimaldict(container[1]), Dict("backend"=>backend(container[1])))
    else
        @test z["inputs"][1]["container"] == merge(minimaldict(container), Dict("backend"=>backend(container)))
    end
    @test z["flow"]["processes"][1]["process"] == "myprocess3"
    @test z["flow"]["processes"][1]["parameters"] == Dict("six"=>6, "seven"=>7)

    rm(io)
end

@testset "default geometry" for cloud in clouds
    container = mkcontainer(cloud, "test-$(uuid4())-cs")
    io = csopen_robust(container, "w", axis_lengths=[10,12,20])
    close(io)

    io = csopen(container)
    g = geometry(io)

    @test g.un == 10
    @test g.vn == 12
    @test g.wn == 20

    rm(io)
end

@testset "iscsdata" for cloud in clouds
    container = mkcontainer(cloud, "test-$(uuid4())-cs")
    io = csopen_robust(container, "w", axis_lengths=[10,12,20])
    close(io)

    @test iscsdata(container)
    rm(io)

    container = mkcontainer(cloud, "test-$(uuid4())-cs")
    @test !iscsdata(container)

    if isa(container, Vector)
        mkpath.(container)
    else
        mkpath(container)
    end
    @test !iscsdata(container)

    if isa(container, Vector)
        write(container[1], "description.json", """{"foo": "bar"}""")
    else
        write(container, "description.json", """{"foo": "bar"}""")
    end
    @test !iscsdata(container)

    if isa(container, Vector)
        rm.(container)
    else
        rm(container)
    end
end