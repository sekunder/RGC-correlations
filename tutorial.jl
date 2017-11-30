helpful_string=""" Triple quotes mean "Read everything between the triple quotes
as a string, verbatim; no escape characters necessary." This is how we'll make
docstrings later. Notice that the closing quotes are on a new line, so there
will be an empty line when this string is printed.
"""

""" Here I define a few variables. Notice that semicolons are only necessary to
separate commands. """

a = 3.0
b = 2; c = 2/3; d = 2//3 # rational numbers wowzers!
x,y = 2.7, -8 # multiple assignments simultaneously!
z = 3.0 + 4im # complex numbers wowzers!

""" Each of these variables has a type, which you can view with `typeof` """
typeof(a)
typeof(b)
typeof(c) # Some languages would return 0, because Int/Int returns Int. Not so in Julia
typeof(d)
typeof(z)
""" We'll get back to those curly braces shortly. """

""" Each of these types is a subtype of the abstract type Number, which has
further subtypes. See https://docs.julialang.org/en/stable/manual/types/ for
more.

Two quick comments, though:
1. There is a "top" abstract type, `Any`. `isa(x, Any)` always returns true.
2. There is a "bottom" type, `Union{}`. `isa(x, Union{})` always returns false.
"""

isa(a, Number) # check if the type of a is a subtype of Number
typeof(a) <: Number # another way to do that
isa(a, Number) == (typeof(a) <: Number) # oh good, these comparisons match
isa(a, Real)
isa(a, Integer)
isa(a, Rational) # Remember, `Rational` and `Integer` are a data types, not the mathematical sets
isa(a, Float64)
isa(a, Float32)

""" Julia does "multiple dispatch", meaning you can define the same function to
behave differently for different input types. Be careful, this is a blessing and
a curse."""
function f(x::Integer)
    return x^2
end
f(x::Complex) = norm(x) # defining f this way is equivalent to the above. Also, the "return" is technically unnecessary.
f(b) # b is an integer, so we get b^2
f(z) # z is complex, so we get the norm
f(a) # But f is not defined for floats, so this gives an error.

""" Suppose we want f to square real numbers and return the norm of complex
numbers. We have three options:

1. f(x::Real) = x^2 and f(x::Complex) = norm(x)
2. f(x::Number) = x^2 and f(x::Complex) = norm(x)
3. f(x) = x^2 and f(x::Complex) = norm(x)

Option 1 is intuitive, but options 2 and 3 work because Julia will find the most
specific version of the function to apply. Try defining f as in option 3 below,
and then see what happens when enter `f(helpful_string)` in the console. """

# Give two definitions of f here!


""" Only abstract types are allowed to have subtypes. Concrete types do not have
subtypes, except in one sense: Parameterization."""

A = [1 2 3; 4 5 6]
B = [1.0 2.0 3.0; 4.0 5.0 6.0]

""" A and B are both `Array`s, so each of these returns true: """
isa(A, Array)
typeof(B) <: Array
""" But this is false, even though Int64 and Float64 are both subtypes of Any and Number: """
typeof(B) <: Array{Any,2}
typeof(B) <: Array{Number}

A == B # should this return true or false? Why?

""" A Julia convention is that functions with . in the name act elementwise, and
function names that end with ! will modify their arguments in-place (as opposed
to returning a brand new object, as is usually the case) """

A .== B # dots indicate that the function operates elementwise
exp.(A)

v = [1,2,3,4,5]
v2 = v # without going into technical details, this is not copying v. It's just creating a new pointer to the array.
w = shuffle(v)
w # this is all shuffled up, but...
v # v is still in order. shuffle returns a *new* vector.
shuffle!(v)
v # but shuffle! modifies it in place, so now v is changed... forever.
w # notice that w is not affected by this.
v2 # but this is, because it's pointing to the same array, which just got changed.

""" Let's dive into how I've used "Unix design philosophy" and Julia's type
system to make my life easier. """
