module BloscExt

using CloudSeis, Blosc

import CloudSeis:BloscCompressor,CACHE_ALL,CACHE_ALL_LEFT_JUSTIFY,CSeis,cache!,cache_from_file!

function CloudSeis.cache_from_file!(io::CSeis{T,N,BloscCompressor}, extentindex, regularize) where {T,N}
    regularize || error("regularize=false is not supported for 'BloscCompressor'")
    @debug "reading and decompressing extent $extentindex..."
    t_read = @elapsed begin
        cdata = read!(io.extents[extentindex].container, io.extents[extentindex].name, Vector{UInt8}(undef, filesize(io.extents[extentindex].container, io.extents[extentindex].name)))
    end

    t_decompress = @elapsed begin
        io_cdata = IOBuffer(cdata; read=true, write=false)
        nbuffers = read(io_cdata, Int)
        buffersize = read!(io_cdata, Vector{Int}(undef, nbuffers))

        io.cache.data = UInt8[]
        for ibuffer = 1:nbuffers
            cbuffer = read!(io_cdata, Vector{UInt8}(undef, buffersize[ibuffer]))
            io.cache.data = [io.cache.data;decompress(UInt8, cbuffer)]
        end
    end
    mb = length(io.cache.data) / 1_000_000
    mb_compressed = length(cdata) / 1_000_000
    mbps_read = mb_compressed/t_read
    mbps_decompress = mb / t_decompress
    mbps = mb / (t_read + t_decompress)
    @debug "...data read and decompressed (effective: $mbps MB/s; decompression: $mbps_decompress MB/s; read: $mbps_read MB/s -- $mb MB; $mb_compressed compressed MB)"
    nothing
end


function CloudSeis.cache_foldmap!(io::CSeis{T,N,BloscCompressor}, extentindex::Integer, force=false) where {T,N}
    if extentindex == io.cache.extentindex && io.cache.type ∈ (CACHE_ALL,CACHE_ALL_LEFT_JUSTIFY) && !force
        return extentindex
    end
    # TODO - we can't to a partial cache here because of how the compression works.
    cache!(io, extentindex, true, force)
    extentindex
end

function Base.flush(io::CSeis{T,N,BloscCompressor}) where {T,N}
    if io.cache.extentindex == 0
        return nothing
    end

    function flush_blosc_buffer_byterange(ibuffer, buffer_size, buffer_remainder)
        isremainder = ibuffer <= buffer_remainder
        firstbyte = (ibuffer - 1)*buffer_size + (isremainder ? ibuffer : buffer_remainder + 1)
        lastbyte = firstbyte + (isremainder ? buffer_size : buffer_size - 1)
        firstbyte,lastbyte
    end

    @debug "compressing and writing extent $(io.cache.extentindex)..."
    t_compress = @elapsed begin
        Blosc.set_num_threads(Sys.CPU_THREADS)
        Blosc.set_compressor(io.cache.compressor.algorithm)

        maxbuffersize = 2_000_000_000 # This is due to a limitation of the Blosc library.
        cachesize = length(io.cache.data)
        nbuffers,remaining_bytes = divrem(cachesize, maxbuffersize)
        remaining_bytes > 0 && (nbuffers += 1)
        buffersize,buffer_remainder = divrem(cachesize, nbuffers)

        cdata = zeros(UInt8, 8*(1 + nbuffers)) # store number of buffers and length of each compressed buffer
        cdata_io = IOBuffer(cdata; read=true, write=true)
        write(cdata_io, nbuffers)
        close(cdata_io)
        for ibuffer = 1:nbuffers
            firstbyte,lastbyte = flush_blosc_buffer_byterange(ibuffer, buffersize, buffer_remainder)
            _data = unsafe_wrap(Array, pointer(io.cache.data)+(firstbyte-1), (lastbyte-firstbyte+1,); own=false)
            _cdata = compress(_data; level=io.cache.compressor.level, shuffle=io.cache.compressor.shuffle)
            cdata = [cdata;_cdata]
            cdata_io = IOBuffer(cdata; read=true, write=true)
            seek(cdata_io, 8*ibuffer)
            write(cdata_io, length(_cdata))
            close(cdata_io)
        end
        close(cdata_io)
    end
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