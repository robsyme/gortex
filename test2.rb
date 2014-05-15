def complement(s)
  s.strip.chars.map do |c|
    case c
    when /[aA]/
      "T"
    when /[cC]/
      "G"
    when /[gG]/
      "C"
    when /[tT]/
      "A"
    else
      c
    end
  end.join
end

def reverse_complement(s)
  complement(s).reverse
end

while ARGF.gets
  puts $_ < reverse_complement($_)
end
