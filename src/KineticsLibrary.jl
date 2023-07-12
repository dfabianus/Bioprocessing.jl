module KineticsLibrary
export monod, blackman
monod(c, q_max, k) = q_max * c / (c + k)
blackman(c, q_max, k) = q_max * minimum([1,c*k])

end