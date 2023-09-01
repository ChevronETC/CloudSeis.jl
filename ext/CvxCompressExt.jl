module CvxCompressExt

using CloudSeis, CvxCompress, ZfpCompression

import CloudSeis:CloudSeisCvxCompressor,CSeis,ZfpCompressor,cachesize,unsafe_foldmap,unsafe_gethdrs,unsafe_gettrcs

function CloudSeis.cache_from_file!(io::CloudSeis.CSeis{T,N,<:CloudSeisCvxCompressor}, extentindex, regularize) where {T,N}
    regularize || error("regularize=false is not supported for 'CloudSeisCvxCompressor'")
    @debug "reading and decompressing extent $extentindex..."
    t_read = @elapsed begin
        cdata = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{UInt8}(undef, filesize(io.extents[extentindex].container, io.extents[extentindex].name)))
    end

    t_decompress = @elapsed begin
        io_cdata = IOBuffer(cdata; read=true, write=false)
        io.cache.data = zeros(UInt8, cachesize(io, extentindex))

        nfmap = read(io_cdata, Int)
        cfmap = read!(io_cdata, Vector{UInt8}(undef, nfmap))
        cfmap = zfp_decompress!(unsafe_foldmap(io, extentindex), cfmap)

        nhdrs = read(io_cdata, Int)
        chdrs = read!(io_cdata, Vector{UInt8}(undef, nhdrs))
        zfp_decompress!(unsafe_gethdrs(io, extentindex), chdrs)
        empty!(chdrs)

        ntrcs_bytes = read(io_cdata, Int)
        ntrcs = div(ntrcs_bytes, 4)
        ctrcs = read!(io_cdata, Vector{UInt32}(undef, ntrcs))
        trcs = unsafe_gettrcs(io, extentindex)

        b3 = min(size(trcs,3), io.cache.compressor.b3)
        if b3 > 1
            b3 = clamp(nextpow(2, b3), 8, 256)
        end

        c = CvxCompress.CvxCompressor((io.cache.compressor.b1,io.cache.compressor.b2,b3), io.cache.compressor.scale)
        CvxCompress.decompress!(trcs, c, ctrcs, ntrcs_bytes)
        empty!(ctrcs)
    end
    mb = length(io.cache.data) / 1_000_000
    mb_compressed = length(cdata) / 1_000_000
    empty!(cdata)
    mbps_read = mb_compressed/t_read
    mbps_decompress = mb / t_decompress
    mbps = mb / (t_read + t_decompress)
    @debug "...data read and decompressed (effective: $mbps MB/s; decompression: $mbps_decompress MB/s; read: $mbps_read MB/s -- $mb MB; $mb_compressed compressed MB)"
    nothing
end

function Base.flush(io::CSeis{T,N,<:CloudSeisCvxCompressor}) where {T,N}
    if io.cache.extentindex == 0
        return nothing
    end

    @debug "compressing and writing extent $(io.cache.extentindex)..."
    t_compress = @elapsed begin
        fmap = unsafe_foldmap(io)
        hdrs = unsafe_gethdrs(io)
        trcs = unsafe_gettrcs(io)

        map(i->begin size(hdrs,i) > typemax(Int32) && error("cvx: each hdrs dimension must be less than $(typemax(Int32)), size(hdrs,$i)=$(size(hdrs,i))") end, 1:ndims(hdrs))

        io_cdata = IOBuffer(;read=false, write=true)

        cfmap = zfp_compress(fmap; nthreads=Sys.CPU_THREADS, write_header=false)
        write(io_cdata, length(cfmap))
        write(io_cdata, cfmap)
        empty!(cfmap)

        chdrs = zfp_compress(hdrs; nthreads=Sys.CPU_THREADS, write_header=false)
        write(io_cdata, length(chdrs))
        write(io_cdata, chdrs)
        empty!(chdrs)

        b3 = min(size(trcs,3), io.cache.compressor.b3)
        if b3 > 1
            b3 = clamp(nextpow(2, b3), 8, 256)
        end

        c = CvxCompress.CvxCompressor((io.cache.compressor.b1,io.cache.compressor.b2,b3),io.cache.compressor.scale)
        ctrcs = Vector{UInt32}(undef, ceil(Int,4*length(trcs)*sizeof(io.traceformat)/4))
        n = CvxCompress.compress!(ctrcs, c, trcs)
        write(io_cdata, n)
        write(io_cdata, view(ctrcs,1:div(n,4)))
        empty!(ctrcs)
    end

    cdata = take!(io_cdata)
    t_write = @elapsed write(io.extents[io.cache.extentindex].container, io.extents[io.cache.extentindex].name, cdata)
    mb = length(io.cache.data)/1_000_000
    mb_compressed = length(cdata)/1_000_000
    mbps_write = mb_compressed/t_write
    mbps_compress = mb/t_compress
    mbps = mb/(t_write+t_compress)
    @debug "...data compressed and wrote (effective: $mbps MB/s ; compression: $mbps_compress MB/s ; write: $mbps_write MB/s -- $mb MB; $mb_compressed compressed MB)"
    nothing
end

end