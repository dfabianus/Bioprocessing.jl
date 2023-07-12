module Bioprocessing

using DataFrames
export HelloWorld, HelloWorld2, HelloWorld3
# Write your package code here.
HelloWorld() = "Hello World!"
HelloWorld2() = "Hello World!2"
HelloWorld3() = "Hello World!3"
#include("observers.jl")
end
