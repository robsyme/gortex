require 'bio'
require 'set'

s = "TACGACGTCGACT"
seq = Bio::Sequence::NA.new(s)

def shear(seq)
  fwd = ("$$$" + seq + "$").
    chars.
    each_cons(4).
    to_set
  rev = ("$$$" + seq.complement + "$").
    chars.
    each_cons(4).
    to_set
  (fwd + rev).
    to_set
    map{|edge| [edge[0..edge.count-2], edge[-1]]}.
    sort_by{|node, edge| node.reverse}
end

def calcL(arr)
  out = []
  arr.each_cons(2) do |first, second|
    out << (first[0] == second[0] ? 0 : 1)
  end
  out << 1
end

def calcF(arr)
  
end


a = shear(seq)
puts a.map{|edge, node| [edge.join, node].join("\t")}
puts
a.zip(calcL(a)) do |both, l|
  puts [l, both[0].join, both[1]].join("\t")
end
