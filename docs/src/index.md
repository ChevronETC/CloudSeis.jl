# CloudSeis.jl

CloudSeis.jl is a Julia library for reading and writing CloudSeis files.  The
CloudSeis data format is designed to be similar to the JavaSeis[1] data format
while adapting to cloud storage (e.g. Azure Blob Storage).  CloudSeis.jl
works-a-round cloud storage latency issues using a caching layer.

[1] https://github.com/ChevronETC/TeaSeis.jl

## Quick start example
```julia
# load the library
using AzStorage, CloudSeis, FolderStorage

# create a new CloudSeis file from an Azure container.  Use a 3D framework (128 samples per trace, 32 traces per frame, and 16 frames per volume)
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container, "w", axis_lengths=[128, 32, 16])

# alternatively, create a new CloudSeis file for POSIX storage
#=
container = Folder("filename-cs")
io = csopen(container, "w", axis_lengths=[128, 32, 16])
=#

# allocate traces and headers for a single frame
trcs, hdrs = allocframe(io)

# populate trcs and hdrs with values
for i = 1:size(io,2)
  set!(prop(io, stockprop[:TRC_TYPE]), hdrs, i, tracetype[:live])
  set!(prop(io, stockprop[:TRACE]), hdrs, i, i)
  set!(prop(io, stockprop[:FRAME]), hdrs, i, 1)
end
rand!(trcs)

# write trcs,hdrs the data-set
writeframe(io, trcs, hdrs)

# close the file (this will also flush buffers to block storage as needed)
close(io)

# open a CloudSeis dataset from an existing container.
io = csopen(container)

# read the first frame:
trcs, hdrs = readframe(io, 1) # out-of-place read
readframe!(io, trcs, hdrs) # in-place read

# access values stored in a trace property in the first trace of the frame
get(prop(io, stockprop[:TRACE]), hdrs, 1)

# close the file
close(io)
```

## csopen / cscreate
A CloudSeis dataset is created/opened with the `csopen` or `cscreate` methods which
return a `CSeis`.  A CloudSeis dataset must have a minimum of 3 dimensions.  For example:
```julia
using AzStorage, CloudSeis

# create a 3D CloudSeis dataset with 10  samples per trace, 11 traces per frame and 12 frames per volume
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
io = csopen(container, "w", axis_lengths=[10,11,12])

# open an existing dataset in read-only model
io = csopen(container, "r")
io = csopen(container) # equivalent to previous line

# open an existing dataset for reading and writing
io = csopen(container, "r+")

# close an open dataset
close(io)

# create a dataset without returning a handle or opening the data
cscreate(container, axis_lengths=[10,11,12])
```
The `cscreate` method is useful, for example, when you need to create the
dataset on the master process and write to it on worker processes.

It is also possible to pass a list of containers to `csopen` and `cscreate`.  In
this case the extents are distributed across all containers, and the meta-information
is in the first container in the list.  This is, typically, used to manually shard data
across multiple Azure storage accounts in order to improve through-put.  For example:
```julia
containers = [AzContainer("filename-cs"; storageaccount="mystorageaccount$i") for i = 1:10]
io = csopen(containers, "w", axis_lengths=[10,11,2])
```
Note that when opening a data-set that is sharded accross multiple containers in "r" or "r+"
modes only the primary container that contains *description.json* needs to be provided.

`csopen` and `cscreate` take a number of keyword arguments to control behavior.  Please
see the **reference** section in this documentation for more information.

# Read/write methods

CloudSeis is a frame based file format.  For `io::CSeis`, allocate memory for a single frame:

```julia
trcs, hdrs = allocframe(io) # allocate memory for traces and headers for a single frame
trcs = allocframetrcs(io) # allocate memory for traces for a single frame
hdrs = allocframehdrs(io) # allocate memory for headers for a single frame
```

Read a frame. `ifrm::Int`, `ivol::Int`, `ihyp::Int` and `i6::Int` must be consistent
with the CloudSeis data context.
```julia
trcs, hdrs = readframe(io, ifrm) # read from 3D data
trcs, hdrs = readframe(io, ifrm, ivol) # read from 4D data
trcs, hdrs = readframe(io, ifrm, ivol, ihyp) # read from 5D data
trcs, hdrs = readframe(io, ifrm, ivol, ihyp, i6) # read from 6D data
...
```

Read a frame (in-place) using pre-allocated memory:
```julia
readframe!(io, trcs, hdrs, ifrm)                # read from 3D data
readframe!(io, trcs, hdrs, ifrm, ivol)          # read from 4D data
readframe!(io, trcs, hdrs, ifrm, ivol, ihyp)    # read from 5D data
readframe!(io, trcs, hdrs, ifrm, ivol, ihyp, i6) # read from 6D data
...
```

Similar method exist for reading only headers:
```julia
hdrs = readframehdrs(io, ifrm) # read from 3D data
hdrs = readframehdrs(io, ifrm, ivol) # read from 4D data
...
readframehdrs!(io, hdrs, ifrm) # in-place read from 3D data
readframehdsr!(io, hdrs, ifrm, ivol) # in-place read from 4D data
...
```
or only traces:
```julia
trcs = readframetrcs(io, ifrm) # read from 3D data
trcs = readframetrcs(io, ifrm, ivol) # read from 3D data
...
readframetrcs!(io, trcs, ifrm) # in-place read from 3D data
readframetrcs!(io, trcs, ifrm, ivol) # in-place read from 4D dadta
...
```

Write a frame.  The frame, volume, etc. indices are determined from the trace properties.
```julia
writeframe(io, trcs, hdrs)
```

To loop over all frames in a dataset of arbitrary dimensions, use
`LogicalIndices`:
```julia
for idx in LogicalIndices(io)
  trcs, hdrs = readframe(io, idx)
end
```
Of course, this can also be used with `readframe!`, `readframetrcs`, `readframetrcs!`,
`readframehdrs` and `readframehdrs!`.

# Fold

Methods for finding the fold of a frame
```julia
fold(io, hdrs) # get fold by examining `hdrs` from a frame
fold(io, ifrm) # get fold for a 3D data-set using the `TraceMap`
fold(io, ifrm, ivol) # get fold for a 4D data-set using the `TraceMap`
...
```

# Alternative read/write methods (N-dimensional slices)

We supply convenience methods for reading and writing arbitrary patches of data.

**Reading**

```julia
trcs,hdrs = read(io, 1:10, 2:3, 4) # read from 3D dataset (frame 4, traces 2-3 and time samples 1-10)
trcs,hdrs = read(io, 1:10, 2:3, 4, :) # read from a 4D data-set (all volumes, frame 4, traces 2-3, and time samples 1-10)
...
read!(io, trcs, hdrs, 1:10, 2:3, 4) # in-place read from 3D data
...
```

Similar methods exist for reading only traces:
```julia
trcs = readtrcs(io, 1:10, 2:3, 4)
readtrcs!(io, trcs, 1:10, 2:3, 4) # in-place version of previous line
```
and only headers:
```julia
hdrs = readhdrs(io, 1:10, 2:3, 4)
readhdrs!(io, hdrs, 1:10, 2:3, 4) # in-place version of previous line
```

**Writing**

```julia
write(io, trcs, hdrs) # trcs::Array{Float32,N}, hdrs::Array{UInt8,N} where N>=3
write(io, trcs, hdrs, 1:10) # same as previous except only time samples 1:10 are written
```

## Alternative write methods for full frames

The first set of APIs are for writing one frame at a time:
```julia
writeframe(io, trcs, ifrm) # write to 3D data
writeframe(io, trcs, ifrm, ivol) # write to 4D data
...
```

The second set of APIs are for writing arbitrary N-dimensional data:
```julia
write(io, trcs, :, 1:10, 3:2:5) # write to 3D data, all samples; traces 1-10; frames 3,5
write(io, trcs, :, 1:10, 3:2:5, 6) # write to 4D data, all samples; traces 1-10; frames 3,5; volume 6
...
```

Please note that in these forms, the `writeframe` and `write` methods will create
headers for you, and populate the `:TRC_TYPE` property along with the properties
corresponding to the trace and frame axes of your data.  In the case of 4D data,
the volume property is also populated, and in the case of 5D data, the volume and
hypercube properties are also populated.

In addition, please note that in the `write` method, `trcs` must have the same
number of dimensions as `io`.  In practice this can be accomplished using `reshape`.
For example if `size(io)=(10,20,3)` and `size(trcs)=(10,)`, then to write `trcs`
to the first trace of the first frame of `io` one could write:

```julia
write(io, rehsape(trcs, 10, 1, 1), :, 1, 1)
```

# Reduction (aggregation) of a CloudSeis data-set
There are scenarios where a CloudSeis data-set will consist of many small
extents.  For example, if one is writing a CloudSeis data-set from many parallel
processes and where each process is responsible for writing a small amount of
data.  Subsequently, one may want to read the data-set from a single instance.
In order to avoid latencies in this subsequent read, CloudSeis.jl provides a
`reduce` method that aggregates frames into fewer extents.  For example:
```julia
using Distributed, AzManagers
addprocs("cbox16", 16)
@everywhere using CloudSeis

# create a data-set with 16 frames and 1 frame per extent
container = AzContainer("mydataset-cs"; storageaccount="mystorageaccount")
cscreate(container, "w", axis_lengths=[100,100,500], extents=[i:i for i=1:16])

# each process writes to its own extent
@sync for (ipid,pid) in enumerate(workers())
  @spawnat pid begin
    io = csopen(Bucket("test-cs"), "r+")
    writeframe(io, rand(Float32,100,100), ipid)
  end
end

# reduce the number of extents such that each extent is about 1GB (GCP only)
reduce(container, mbytes_per_extent=1000)

rmprocs(workers())
```
Please see the **reference** section of this documentation for more information.

# Trace Properties

The CloudSeis data format does not specify any trace properties.  However, there
are commonly used (**stock**) properties (see `src/stockprops.jl`).  It is
unusual when a stock property does not suit your needs.  But, if need be, you
can define a custom property using the `TracePropertyDef` constructor:

```julia
pdef = TracePropertyDef("label", "description", Float32)
pdef = TracePropertyDef("label", "description", Vector{Float32}, 2)
```

The arguments to `TracePropertyDef` are the `label`, `description`, `type`, and,
optionally, the **number of elements** stored in the property. The stock
properties are defined in `src/stockprops.jl` using a Julia
dictionary: `stockprop`.  For example, access a stock definition for the `TRACE`
property:
```julia
pdef = stockprop[:TRACE]
```

Given a CloudSeis file `io::CSeis` and a stock definition, we can access the
corresponding property of a CloudSeis file:

```julia
p = prop(io, pdef) # access using a `TracePropertyDef`
p = prop(io, "TRACE") # alternatively, access using the trace property definition label
p = prop(io, "TRACE", Int32) # type-stable version of previous line
```

Given, additionally, a frame of headers `hdrs::Array{UInt8,2}`, we can get and
set the values stored in a property:

```julia
@show get(p, hdrs, 1) # get trace property value for the first traces in `hdrs`
set!(p, hdrs, 1, 5) # set the first header in `hdrs` to 5
writeframe(io, trcs, hdrs) # the CloudSeis file does not know about the updated header until you call `writeframe`
```
In the above code listing `trcs` is of type `Array{Float32,2}`.

## TRC_TYPE

The `TRC_TYPE` property is used to indicate if a trace is dead, live or auxiliary
within any given frame.  It is stored as an `Int`.  We provide a second
dictionary to map between the `Int` and human readable code:

```julia
tracetype[:live]
tracetype[:dead]
tracetype[:aux]
```

For example,

```julia
io = csopen("file-cs", "r")
trcs, hdrs = readframe(io, 1)
prop_trctype = prop(io, stockprop[:TRC_TYPE])
for i=1:size(hdrs,2)
    if get(prop_trctype, hdrs, i) == tracetype[:live]
        write(STDOUT, "trace $(i) is a live trace\n")
    elseif get(prop_trctype, hdrs, i) == tracetype[:dead]
        write(STDOUT, "trace $(i) is a dead trace\n")
    elseif get(prop_trctype, hdrs, i) == tracetype[:aux]
        write(STDOUT, "trace $(i) is a aux trace\n")
    end
end
close(io)
```

# Data properties

CloudSeis.jl provides support for storing custom data properties.  This is
accomplished by passing an array of `DataProperty`'s to the `csopen` function.
For example, a data property could be defined as:


```julia
p = DataProperty("Survey date", 120977)
```

# Geometry

CloudSeis.jl provides support for storing survey geometry using three-points to
define rotated/translated coordinate system.

```julia
geom = Geometry(u1=1,un=2,v1=1,vn=2,w1=1,wn=2,ux=1.0,uy=0.0,uz=0.0,vx=0.0,vy=1.0,vz=0.0,wx=0.0,wy=0.0,wz=1.0)
```

where `(ox,oy,oz)` is the origin, `(ux,uy,uz)` is a vector to define the end of
the `u-axis` (e.g. cross-line axis), `(vx,vy,vz)` is the end of the `v-axis`
(e.g. the in-line axis), and `(wx,wy,wz)` is the end of the `w-axis`
(e.g. the depth axis).  `(u1,un)` are the first and last bin indices along the
`u-axis`, `(v1,vn)` are the first and last bin indices along the `v-axis`, and
`(w1,wn)` are the first and last bin indices along the `w-axis`.  CloudSeis.jl
does not provide any tools for using this geometry to manipulate trace coordinates.
I would recommend that this functionality be put into a separate package.

# History

CloudSeis.jl provides support for storing processing history by recording the input data-set(s) and steps in the processing flow that resulted in the data-set.  A step is defined as the process (program) that was run as well as the input parameters for that process.  The history is recursive in the sense that the history of input data-sets are embedded.

Example of creating an new history dictionary:

```julia
h = history!(flow_parameters=Dict("one"=>1))
h = history!(h, process="myprocess1", process_parameters=Dict("two"=>2,"three"=>3))
h = history!(h, process="myprocess2", process_parameters=Dict("four"=>4,"five"=>5))
```

We can then use that history in the construction of a new CloudSeis data-set,
```julia
io = csopen(Folder("file-1-cs"), "w", axis_lengths=[12,11,10], history=h)
history(io) # outputs the history as a dictionary
close(io)
```

Finally, we can use embed this history of `file-1-cs` into a new data-set,
```julia
io = csopen(Folder("file-1-cs"))
h = history!(process="myprocess3", process_parameters=Dict("eight"=>8,"nine"=>9), histories = [Folder("file-1-cs")])
close(io)
io = csopen(Folder("file-2-cs"), "w"; axis_lengths=[12,11,10], history=h)
using JSON
print(json(history(io), 1))
```
Then the history structure is:
```json
{
 "flow": {
  "parameters": {},
  "processes": [
   {
    "parameters": {
     "eight": 8,
     "nine": 9
    },
    "process": "myprocess3"
   }
  ]
 },
 "inputs": [
  {
   "history": {
    "flow": {
     "parameters": {
      "one": 1
     },
     "processes": [
      {
       "parameters": {
        "two": 2,
        "three": 3
       },
       "process": "myprocess1"
      },
      {
       "parameters": {
        "four": 4,
        "five": 5
       },
       "process": "myprocess2"
      }
     ]
    }
   },
   "container": {
    "foldername": "/home/tqff/.julia/dev/CloudSeis/file-1-cs"
   }
  }
 ]
}
```

# Convenience methods and dictionaries

For convenience and consistency, we supply several dictionaries.  In addition to
the dictionary for trace property definitions and trace type (both described above),
there are dictionaries for **data domain** `stockdomain`, **units**
`stockunit`, and **data type** `stockdatatype`.  All of these are listed in
`src\stockprops.jl`.


Example usage within the csopen method:


```julia
io = csopen(Bucket("file-cs"), "w", axis_lengths=[12,11,10], axis_units=[stockunit[:SECONDS], stockunit[:METERS], stockunit[:METERS]], axis_domains=[stockdomain[:TIME], stockdomain[:SPACE], stockdomain[:SPACE], datatype=stockdatatype[:SOURCE])
```

Several convenience methods are supplied for querying `io::CSeis`:


```julia
ndims(io)              # returns `Int`, number of dimensions in the CloudSeis dataset
length(io)             # returns `Int`, the number of frames in the CloudSeis dataset, equivalent to `prod(size(io)[3:end])`
size(io)               # returns `NTuple{Int}`, size of CloudSeis dataset
size(io,i)             # returns `Int`, size of CloudSeis dataset along dimension `i::Int`
props(io)              # returns `NTuple{TraceProperty}`, trace property along all dimensions
props(io,i)            # returns `TraceProperty`, trace property along dimension `i::Int`
propdefs(io)           # returns `NTuple{TracePropertyDef}`, trace property definition along all dimensions
propdefs(io,i)         # returns `TracePropertyDef`, trace property along dimension `i::Int`
units(io)              # returns `NTuple{String}`, units along all dimensions
units(io,i)            # returns `String, unit along dimension `i::Int`
domains(io)            # returns `NTuple{String}`, data domains along all dimensions
domains(io,i)          # returns `String`, data domain along dimension `i::Int`
pstarts(io)            # returns `NTuple{Float64}`, physical starts along all dimensions
pstarts(io,i)          # returns `Float64`, physical start along dimension `i::Int`
pincs(io)              # returns `NTuple{Float64}`, physical increments along all dimensions
pincs(io,i)            # returns `Float64`, physical increment along dimension `i::Int`
lstarts(io)            # returns `NTuple{Int}`, logical starts along all dimensions
lstarts(io,i)          # returns `Int`, logical start along dimension `i::Int`
lincs(io)              # returns `NTuple{Int}`, logical increments along all dimensions
lincs(io,i)            # returns `Int`, logical increment along dimension `i::Int`
in(prop,io)            # returns true if the trace property `prop` exists in `io` --  `prop` can be of types `::TraceProperty`, `::TracePropertyDef`, or `::String`
dataproperty(io,nm)    # returns the value held in the data property: `nm::String`
hasdataproperty(io,nm) # returns true if the data property corresponding to label `nm::String` is in `io::CSeis`
geometry(io)           # returns `Geometry`, the stored geometry of the dataset.  If no geometry is stored, `nothing` is returned
```

Convenience methods are supplied for manipulating `io::CSeis`:

```julia
rm(io) # remove (delete) the file and all of its extent files and secondary folders
empty!(io)  # remove extends and secondary folders, but keep meta-data
cp(src, dst) # create a new CloudSeis file `dst::AbstractString` that is a copy of `src::CSeis`, optional named argument: `secondaries=` - change file extents location
mv(src, dst)  # move a CloudSeis file to `dst::AbstractString` from `src::CSeis`, optional named argument: `secondaries=` - change file extents location
copy!(io, hdrs, io1, hdrs1) # (not-implemented) copy values from `hdrs1::Array{UInt8,2}` to `hdrs::Array{UInt8,2}`
reduce(io::CSeis[; mbytes_per_extent=1000, frames_per_extent=0, maxinstances=100, instancetemplate="jbox08", instancegroup="jbox", retries=0, maxerrors=Inf]) # reduce the number of extents used to store data, GCP only
```
