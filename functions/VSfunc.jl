
function getAMwave(clickfreq=2,sr=1000)

    t  = 0:1/sr:1 # time course
    wave = sin.((2*π*clickfreq*t).+3π/2)
    return wave

end