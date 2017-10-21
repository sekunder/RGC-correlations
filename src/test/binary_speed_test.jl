
# Quick check seems to show that, unsurprisingly, the chunks method for
# bitvectors is faster than anything. For Vector{Bool}, it seems that the dot
# product method is faster to a point, but the speed gained from the bitshifting
# method is not really all that much for higher numbers of bits. So, dot product
# is probably the way to go.
binary_to_int_bv(x::BitVector) = Int(x.chunks[1])
binary_to_int_bs(x) = sum([x[i] << (i-1) for i = 1:length(x)])
binary_to_int_dot(x) = dot([2^i for i = 0:(length(x) - 1)], x)

n_trials = 100000; n_bits = [10, 20, 30, 40]

for k in n_bits
  println("Testing k = $k bits. Generating $n_trials random binary vectors, both as Vector{Bool} and as BitVector")
  bool_vecs = [rand(Bool, k) for t = 1:n_trials]
  bit_vecs = map(BitVector, bool_vecs)
  println("\tTesting bitshifting method, Vector{Bool}")
  @time bool_bs = [binary_to_int_bs(x) for x in bool_vecs]
  println()

  # println("\tTesting bitshifting method, BitVector")
  # @time bit_bs = [binary_to_int_bs(x) for x in bit_vecs]
  # println()

  println()


  println("\tTesting dot product method, Vector{Bool}")
  @time bool_dot = [binary_to_int_dot(x) for x in bool_vecs]
  println()

  # println("\tTesting dot product method, BitVector")
  # @time bit_dot = [binary_to_int_dot(x) for x in bit_vecs]
  # println()

  println()

  # println("\tTesting chunks methd (BitVector only)")
  # @time bit_chunk = [binary_to_int_bv(x) for x in bit_vecs]
  # println()

  println("\tChecking that all values agree")
  println("\tbool_bs / bool_dot : $(all(bool_bs .== bool_dot))")
  # println("\tbool_bs  / bit_bs   : $(all(bool_bs .== bit_bs))")
  # println("\tbit_bs   / bool_dot : $(all(bit_bs .== bool_dot))")
  # println("\tbool_dot / bit_dot  : $(all(bool_dot .== bit_dot))")
  # println("\tbit_dot  / bit_chunk: $(all(bit_dot .== bit_chunk))")
end
