module ZfpExt

using CloudSeis, ZfpCompression

import CloudSeis:CACHE_ALL,CACHE_ALL_LEFT_JUSTIFY,CACHE_FOLDMAP,CSeis,CloudSeisCvxCompressor,ZfpCompressor,cachesize,cachesize_foldmap,kwargs,partialcache_error,unsafe_foldmap,unsafe_gethdrs,unsafe_gettrcs

function CloudSeis.cache_from_file!(io::CSeis{T,N,<:ZfpCompressor}, extentindex, regularize) where {T,N}
    regularize || error("regularize=false is not supported for 'ZfpCompressor'")
    @debug "reading and decompressing extent $extentindex..."
    t_read = @elapsed begin
        cdata = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{UInt8}(undef, filesize(io.extents[extentindex].container, io.extents[extentindex].name)))
    end

    t_decompress = @elapsed begin
        io_cdata = IOBuffer(cdata; read=true, write=false)
        io.cache.data = zeros(UInt8, cachesize(io, extentindex))

        nfmap = read(io_cdata, Int)
        cfmap = read!(io_cdata, Vector{UInt8}(undef, nfmap))
        zfp_decompress!(unsafe_foldmap(io, extentindex), cfmap)
        empty!(cfmap)

        nhdrs = read(io_cdata, Int)
        chdrs = read!(io_cdata, Vector{UInt8}(undef, nhdrs))
        zfp_decompress!(unsafe_gethdrs(io, extentindex), chdrs)
        empty!(chdrs)

        ntrcs = read(io_cdata, Int)
        ctrcs = read!(io_cdata, Vector{UInt8}(undef, ntrcs))
        zfp_decompress!(unsafe_gettrcs(io, extentindex), ctrcs; kwargs(io.cache.compressor)...)
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

function CloudSeis.read_foldmap!(io::Union{CSeis{T,N,<:ZfpCompressor},CSeis{T,N,CloudSeisCvxCompressor}}, extentindex::Integer, fmap::Vector{UInt8}; serial=false) where {T,N}
    try
        nfmap = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{Int}(undef, 1); serial=false)[1]
        cfmap = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{UInt8}(undef, nfmap); offset=sizeof(Int), serial)
        _fmap = unsafe_wrap(Array, convert(Ptr{Int}, pointer(fmap)), length(io.extents[extentindex].frameindices); own=false)
        zfp_decompress!(_fmap, cfmap)
    catch e
        fmap .= 0
        if !isa(e, FileDoesNotExistError)
            throw(e)
        end
    end
    fmap
end

function CloudSeis.cache_foldmap!(io::Union{CSeis{T,N,<:ZfpCompressor},CSeis{T,N,CloudSeisCvxCompressor}}, extentindex::Integer, force=false) where {T,N}
    io.mode == "r" || partialcache_error()

    if extentindex == io.cache.extentindex && io.cache.type ∈ (CACHE_ALL,CACHE_ALL_LEFT_JUSTIFY) && !force
        return extentindex
    end

    if isfile(io.extents[extentindex].container, io.extents[extentindex].name)
        @debug "reading foldmap for extent $extentindex from storage..."
        t_read = @elapsed begin
            nfmap = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{Int}(undef, 1))[1]
            cfmap = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{UInt8}(undef, nfmap); offset=sizeof(Int))
        end
        t_decompress = @elapsed begin
            io.cache.data = zeros(UInt8, cachesize_foldmap(io, extentindex))
            zfp_decompress!(unsafe_foldmap(io, extentindex), cfmap)
        end
        mb = length(io.cache.data) / 1_000_000
        mb_compressed = length(cfmap) / 1_000_000
        empty!(cfmap)
        mbps_read = mb_compressed/t_read
        mbps_decompress = mb / t_decompress
        mbps = mb / (t_read + t_decompress)
        @debug "...foldmap read and decompressed (effective: $mbps MB/s; decompression: $mbps_decompress MB/s; read: $mbps_read MB/s -- $mb MB; $mb_compressed compressed MB)"
    else
        @debug "creating cache..."
        io.cache.data = zeros(UInt8, cachesize_foldmap(io, extentindex))
        @debug "...done creating cache."
    end
    io.cache.extentindex = extentindex
    io.cache.type = CACHE_FOLDMAP

    extentindex
end

function Base.flush(io::CSeis{T,N,<:ZfpCompressor}) where {T,N}
    if io.cache.extentindex == 0
        return nothing
    end

    @debug "compressing and writing extent $(io.cache.extentindex)..."
    t_compress = @elapsed begin
        fmap = unsafe_foldmap(io)
        hdrs = unsafe_gethdrs(io)
        trcs = unsafe_gettrcs(io)

        map(i->begin size(hdrs,i) > typemax(Int32) && error("zfp: each hdrs dimension must be less than $(typemax(Int32)), size(hdrs,$i)=$(size(hdrs,i))") end, 1:ndims(hdrs))
        map(i->begin size(trcs,i) > typemax(Int32) && error("zfp: each trcs dimension must be less than $(typemax(Int32)), size(hdrs,$i)=$(size(trcs,i))") end, 1:ndims(trcs))

        io_cdata = IOBuffer(;read=false, write=true)

        cfmap = zfp_compress(fmap; nthreads=Sys.CPU_THREADS, write_header=false)
        write(io_cdata, length(cfmap))
        write(io_cdata, cfmap)
        empty!(cfmap)

        chdrs = zfp_compress(hdrs; nthreads=Sys.CPU_THREADS, write_header=false)
        write(io_cdata, length(chdrs))
        write(io_cdata, chdrs)
        empty!(chdrs)

        ctrcs = zfp_compress(trcs; nthreads=Sys.CPU_THREADS, write_header=false, kwargs(io.cache.compressor)...)
        write(io_cdata, length(ctrcs))
        write(io_cdata, ctrcs)
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