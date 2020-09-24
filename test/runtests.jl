using AzSessions, AzStorage, CloudSeis, Dates, Distributed, FolderStorage, HTTP, Random, Test

const client_id = get(ENV, "CLIENT_ID", "")
const client_secret = get(ENV, "CLIENT_SECRET", "")
const storageaccount = get(ENV, "STORAGEACCOUNT1", "unittester1")

#const session_az = AzSession(;protocal=AzClientCredentials, client_id=client_id, client_secret=client_secret, resource="https://storage.azure.com")

const storageaccount1 = get(ENV, "STORAGEACCOUNT1", "unittester")
const storageaccount2 = get(ENV, "STORAGEACCOUNT2", "unittester2")
const storageaccounts = [storageaccount1, storageaccount2]

# function new_container_az(foldername)
#     AzContainer(foldername, session=session_az, storageaccount=storageaccounts[1])
# end

# function new_container_az2(foldername)
#     [AzContainer(foldername, session=session_az, storageaccount=sa) for sa in storageaccounts]
# end

function new_container_posix(foldername)
    Folder(foldername)
end

function new_container(cloud, foldername)
    if cloud == "AZ"
        return new_container_az(foldername)
    elseif cloud == "AZ2"
        return new_container_az2(foldername)
    elseif cloud == "POSIX"
        return new_container_posix(foldername)
    end
end

# const clouds = ("AZ","AZ2","POSIX")
const clouds = ("POSIX",)
const compressors = ("none","blosc")

@testset "CloudSeis, cloud=$cloud, compresser=$compressor" for cloud in clouds, compressor in compressors
    @testset "allocframe" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+0)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)
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
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)

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

        io = csopen(new_container(cloud, "test-$r-cs"))
        t,h = readframe(io, 1)
        close(io)
        @test t ≈ x
        rm(io)
    end

    @testset "writeframe/readframe with index" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+2)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], dataproperties=[DataProperty("P",1)], compressor=compressor)

        x = rand(Float32,10,11)

        t,h = allocframe(io)
        t .= x

        writeframe(io, t, 1)
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"))
        t,h = readframe(io, 1)
        close(io)
        @test t ≈ x
        rm(io)
    end

    @testset "write" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+3)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"))
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
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], frames_per_extent=1, compressor=compressor)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"))
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
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"))
        _x = readtrcs(io, :, :, :)
        close(io)
        rm(io)

        @test x ≈ _x
    end

    @testset "readhdrs" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+6)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"))
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
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"))
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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+6)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)

        x = rand(Float32,10,11,12)
        write(io, x, :, :, :)
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"))
        f = [fold(io,i) for i=1:12]
        close(io)
        rm(io)

        @test f == [11 for i=1:12]
    end

    @testset "similarto" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+7)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w",
            axis_pincs=[0.1,0.2,0.3],
            axis_lengths=[10,11,12],
            dataproperties=[DataProperty("P",1)],
            compressor=compressor)
        close(io)

        _io = csopen(new_container(cloud, "test-$r-sim-cs"), "w", similarto=new_container(cloud, "test-$r-cs"))
        @test size(_io) == (10,11,12)
        @test pincs(_io)[1] ≈ 0.1
        @test pincs(_io)[2] ≈ 0.2
        @test pincs(_io)[3] ≈ 0.3

        rm(io)
        rm(_io)
    end

    @testset "similarto, changing extents" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+7)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w",
            axis_pincs=[0.1,0.2,0.3],
            axis_lengths=[10,11,12],
            dataproperties=[DataProperty("P",1)],
            compressor=compressor)
        @test length(io.extents) == 1
        close(io)

        _io = csopen(new_container(cloud, "test-$r-sim-cs"), "w", similarto=new_container(cloud, "test-$r-cs"), frames_per_extent=1)
        @test size(_io) == (10,11,12)
        @test pincs(_io)[1] ≈ 0.1
        @test pincs(_io)[2] ≈ 0.2
        @test pincs(_io)[3] ≈ 0.3

        @test length(_io.extents) == 12 
    end

    @testset "propdefs" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+8)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)
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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+9)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)
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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+10)))
        pdef = TracePropertyDef("X","XX",Vector{Float64},2)
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], tracepropertydefs=[pdef], compressor=compressor)
        t,h = allocframe(io)
        for i = 1:11
            set!(prop(io,"TRACE"), h, i, i)
            set!(prop(io,"FRAME"), h, i, 1)
            set!(prop(io,"TRC_TYPE"), h, i, tracetype[:live])
            set!(prop(io,"X"), h, i, i*[1.0,2.0])
        end
        writeframe(io,t,h)
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"))
        t,h = readframe(io,1)

        for i = 1:11
            @test get(prop(io,"X"), h, i) ≈ i*[1.0,2.0]
        end
        close(io)
        rm(io)
    end

    @testset "pincs" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+11)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_pincs=[1.0,2.0,3.0], compressor=compressor)
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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+12)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_pstarts=[1.0,2.0,3.0], compressor=compressor)
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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+13)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_units=["X","Y","Z"], compressor=compressor)
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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+14)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], axis_domains=["X","Y","Z"], compressor=compressor)
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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+15)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], geometry=g, compressor=compressor)

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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+16)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)

        @test in(stockprop[:TRACE], io) == true
        @test in(stockprop[:CDP],io) == false

        close(io)
        rm(io)
    end

    @testset "dataproperty" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+17)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], dataproperties=[DataProperty("X",1)], compressor=compressor)

        @test dataproperty(io, "X") == 1

        close(io)
        rm(io)
    end

    @testset "hasdataproperty" begin
        sleep(1)
        r = lowercase(randstring(MersenneTwister(millisecond(now())+18)))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], dataproperties=[DataProperty("X",1)], compressor=compressor)

        @test hasdataproperty(io, "X") == true
        @test hasdataproperty(io, "Y") == false

        close(io)
        rm(io)
    end

    @testset "copy!" begin
        sleep(1)
        r1 = lowercase(randstring(MersenneTwister(millisecond(now())+19),4))
        io1 = csopen(new_container(cloud, "test-$r1-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)
        r2 = lowercase(randstring(MersenneTwister(millisecond(now())+20),4))
        io2 = csopen(new_container(cloud, "test-$r2-cs"), "w", axis_lengths=[10,11,12], compressor=compressor)

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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+21),4))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,50], frames_per_extent=2, compressor=compressor)

        @test length(io.extents) == 25

        T = rand(Float32,10,11,50)
        H = rand(UInt8,headerlength(io),11,50)
        for i = 1:50
            writeframe(io, T[:,:,i], i)
            readframehdrs!(io, view(H,:,:,i), i)
        end
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"), compressor=compressor)
        reduce(io)
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"), compressor=compressor)
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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+22),4))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,50], frames_per_extent=2, compressor=compressor)

        @test length(io.extents) == 25

        T = rand(Float32,10,11,50)
        H = rand(UInt8,headerlength(io),11,50)
        for i = 1:50
            writeframe(io, T[:,:,i], i)
            readframehdrs!(io, view(H,:,:,i), i)
        end
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"))
        addprocs(5)
        @everywhere using AzStorage, CloudSeis, FolderStorage
        reduce(io; frames_per_extent=10)
        rmprocs(workers())
        close(io)

        io = csopen(new_container(cloud, "test-$r-cs"))
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
        r = lowercase(randstring(MersenneTwister(millisecond(now())+23),4))
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], force=true, compressor=compressor)
        writeframe(io, rand(Float32,10,11), 1)
        close(io)
        io = csopen(new_container(cloud, "test-$r-cs"), "w", axis_lengths=[10,11,12], force=true, compressor=compressor)
        writeframe(io, rand(Float32,10,11), 1)
        close(io)
    end
end
