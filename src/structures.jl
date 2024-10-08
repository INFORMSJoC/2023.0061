mutable struct fnode
    id::Int
    name :: String
    ntr::Int
end
mutable struct node
    n :: fnode
    t :: Int
    id :: Int
end
mutable struct arc
    tail :: node
    head :: node
    tt :: Int
    id :: Int
end
mutable struct arck
    tail :: node
    head :: node
    k :: Int
    id :: Int
end
mutable struct commk
    ok :: fnode
    tp :: Int
    td :: Int
    id :: Int
end
mutable struct scommk
    ok :: fnode
    list
    tp :: Int
    td :: Int
    tpr :: Int
    tdr :: Int
    id :: Int
end
function Base.push!(fnodes::node,fnode::node)
    push!(fnodes,fnode)
end
# function Base.push!(nodes::node,node::node)
#     push!(nodes,node)
# end
function Base.push!(arcs::arc,arc::arc)
    push!(arcs,arc)
end
function Base.push!(arcks::arck,arck::arck)
    push!(arcks,arck)
end
# function Base.push!(commks::arck,commk::arck)
#     push!(commks,commk)
# end
# function Base.push!(scommks::arck,scommk::arck)
#     push!(scommks,scommk)
# end
mutable struct fnode2
    id:: Int
    name :: String
    p:: Int
    ag :: Int
    cot :: Int64
    bg :: Int64
    delta_g :: Int64
    tmax_g :: Int64
    temp_c :: Int64
    temp_s :: Int64
    ntr :: Int
    TMAX :: Int
end
mutable struct fnode3
    id:: Int
    lat
    long
end
# function Base.push!(fnodes2::node,fnode2::node)
#     push!(fnodes2,fnode2)
# end

mutable struct v_i
    id :: Int64
    cc :: Int64
end