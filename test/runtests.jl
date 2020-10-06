using AzSessions, AzStorage, CloudSeis, Dates, Distributed, FolderStorage, HTTP, JSON, Random, Test

credentials = JSON.parse(ENV["AZURE_CREDENTIALS"])
AzSessions.write_manifest(;client_id=credentials["clientId"], client_secret=credentials["clientSecret"], tenant=credentials["tenantId"])

session = AzSession(;protocal=AzClientCredentials, client_id=credentials["clientId"], client_secret=credentials["clientSecret"], resource="https://storage.azure.com/")

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
           io = csopen(containers, mode; kwargs...)
           break
        catch e
            @warn "caught excption in csopen, sleeping for 60 seconds"
            showerror(stdout, e)
            sleep(60)
            if time() - tic > timeout
                throw(e)
            end
        end
    end
    io
end

const clouds = (Azure, Azure2, POSIX)

@testset "CloudSeis, cloud=$cloud" for cloud in clouds
    @testset "allocframe" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+0)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12])
        t,h = allocframe(io)

        @test size(t) == (10,11)
        @test eltype(t) == Float32

        @test size(h) == (io.hdrlength,11)
        @test eltype(h) == UInt8

        close(io)
        rm(io)
    end

    @testset "writeframe/readframe with headers" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+1)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12])

        x = rand(Float32,10,11)

        t,h = allocframe(io)
        t .= x

        for i = 1:size(h,2)
            set!(prop(io,io.axis_propdefs[2]), h, i, i)
            set!(prop(io,io.axis_propdefs[3]), h, i, 1)
            set!(prop(io,stockprop[:TRC_TYPE]), h, i, tracetype[:live])
        end
        writeframe(io,t,h)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        t,h = readframe(io, 1)
        close(io)
        @test t ≈ x
        rm(io)
    end

    @testset "writeframe/readframe with index" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+2)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], dataproperties=[DataProperty("P",1)])

        x = rand(Float32,10,11)

        t,h = allocframe(io)
        t .= x

        writeframe(io, t, 1)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        t,h = readframe(io, 1)
        close(io)
        @test t ≈ x
        rm(io)
    end

    @testset "write" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+3)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12])

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        for ifrm = 1:12
            _x = readframetrcs(io, ifrm)
            @test _x ≈ x[:,:,ifrm]
        end
        close(io)
        rm(io)
    end

    @testset "write, multiple extents" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+4)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], frames_per_extent=1)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        for ifrm = 1:12
            _x = readframetrcs(io, ifrm)
            @test _x ≈ x[:,:,ifrm]
        end
        close(io)
        rm(io)
    end

    @testset "readtrcs" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+5)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12])

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        _x = readtrcs(io, :, :, :)
        close(io)
        rm(io)

        @test x ≈ _x
    end

    @testset "readhdrs" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+6)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12])

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        h = readhdrs(io, :, :, :)
        close(io)
        rm(io)

        for ifrm = 1:12, itrc = 1:11
            @test get(prop(io, "TRACE"), view(h, :, :, ifrm), itrc) ≈ itrc
            @test get(prop(io, "FRAME"), view(h, :, :, ifrm), itrc) ≈ ifrm
        end
    end

    @testset "readhdrs!" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+7)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12])

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        h = readhdrs(io, :, :, :)
        close(io)
        rm(io)

        for ifrm = 1:12, itrc = 1:11
            @test get(prop(io, "TRACE"), view(h, :, :, ifrm), itrc) ≈ itrc
            @test get(prop(io, "FRAME"), view(h, :, :, ifrm), itrc) ≈ ifrm
        end
    end

    @testset "partial read for foldmap" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+8)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12])

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        f = [fold(io,i) for i=1:12]
        close(io)
        rm(io)

        @test f == [11 for i=1:12]
    end

    @testset "similarto" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+9)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w",
            axis_pincs=[0.1,0.2,0.3],
            axis_lengths=[10,11,12],
            dataproperties=[DataProperty("P",1)])
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
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+10)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w",
            axis_pincs=[0.1,0.2,0.3],
            axis_lengths=[10,11,12],
            dataproperties=[DataProperty("P",1)])
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
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+11)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12])
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
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+12)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12])
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
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+13)))
        pdef = TracePropertyDef("X","XX",Vector{Float64},2)
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], tracepropertydefs=[pdef])
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

    @testset "pincs" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+14)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_pincs=[1.0,2.0,3.0])
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
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+15)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_pstarts=[1.0,2.0,3.0])
        @test pstarts(io,1) ≈ 1.0
        @test pstarts(io,2) ≈ 2.0
        @test pstarts(io,3) ≈ 3.0
        @test pstarts(io)[1] ≈ 1.0
        @test pstarts(io)[2] ≈ 2.0
        @test pstarts(io)[3] ≈ 3.0
        close(io)
        rm(io)
    end

    @testset "units" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+16)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_units=["X","Y","Z"])
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
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+17)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_domains=["X","Y","Z"])
        @test domains(io,1) == "X"
        @test domains(io,2) == "Y"
        @test domains(io,3) == "Z"
        @test domains(io)[1] == "X"
        @test domains(io)[2] == "Y"
        @test domains(io)[3] == "Z"
        close(io)
        rm(io)
    end

    @testset "geometry" begin
        sleep(1)
        g = Geometry(
            ox=1.0,oy=2.0,oz=3.0,
            ux=4.0,uy=5.0,uz=6.0,
            vx=7.0,vy=8.0,vz=9.0,
            wx=10.0,wy=11.0,wz=12.0,
            u1=1,un=2,
            v1=3,vn=4,
            w1=5,wn=6)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+18)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], geometry=g)

        _g = geometry(io)
        @test g.ox ≈ 1.0
        @test g.oy ≈ 2.0
        @test g.oz ≈ 3.0
        @test g.ux ≈ 4.0
        @test g.uy ≈ 5.0
        @test g.uz ≈ 6.0
        @test g.vx ≈ 7.0
        @test g.vy ≈ 8.0
        @test g.vz ≈ 9.0
        @test g.wx ≈ 10.0
        @test g.wy ≈ 11.0
        @test g.wz ≈ 12.0
        @test g.u1 == 1
        @test g.un == 2
        @test g.v1 == 3
        @test g.vn == 4
        @test g.w1 == 5
        @test g.wn == 6

        close(io)
        rm(io)
    end

    @testset "in" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+19)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12])

        @test in(stockprop[:TRACE], io) == true
        @test in(stockprop[:CDP],io) == false

        close(io)
        rm(io)
    end

    @testset "dataproperty" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+20)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], dataproperties=[DataProperty("X",1)])

        @test dataproperty(io, "X") == 1

        close(io)
        rm(io)
    end

    @testset "hasdataproperty" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+21)))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], dataproperties=[DataProperty("X",1)])

        @test hasdataproperty(io, "X") == true
        @test hasdataproperty(io, "Y") == false

        close(io)
        rm(io)
    end

    @testset "copy!" begin
        sleep(1)
        r1 = lowercase(randstring(MersenneTwister(millisecond(now())+22),4))
        io1 = csopen_robust(mkcontainer(cloud, "test-$r1-cs"), "w", axis_lengths=[10,11,12])
        sleep(1)
        r2 = lowercase(randstring(MersenneTwister(millisecond(now())+23),4))
        io2 = csopen_robust(mkcontainer(cloud, "test-$r2-cs"), "w", axis_lengths=[10,11,12])

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
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+24),4))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,50], frames_per_extent=2)

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
            @test t ≈ T[:,:,i]
            @test h == H[:,:,i]
        end

        @test length(io.extents) == 1

        rm(io)
    end

    @testset "reduce, parallel" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+25),4))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,50], frames_per_extent=2)

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
        @everywhere using AzStorage, CloudSeis, FolderStorage
        reduce(io; frames_per_extent=10)
        rmprocs(workers())
        close(io)

        io = csopen(mkcontainer(cloud, "test-$r-cs"))
        for i = 1:50
            t,h = readframe(io, i)
            @test t ≈ T[:,:,i]
            @test h == H[:,:,i]
        end

        @test length(io.extents) == 5

        rm(io)
    end

    @testset "robust cscreate" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+26),4))
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], force=true)
        writeframe(io, rand(Float32,10,11), 1)
        close(io)
        io = csopen_robust(mkcontainer(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], force=true)
        writeframe(io, rand(Float32,10,11), 1)
        close(io)
    end
end
