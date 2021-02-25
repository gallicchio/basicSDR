# basicSDR

﻿SDR Complex Mixing, Sampling, Fourier, Zero-IF Quadrature Direct-Conversion

1. sampling
    Show cosine
        Slider to higher freq
        looks like crap with fewer than 10 samples per cycle
        but spectrum looks ok even at high freq
        max freq at half of sample rate
        above that, get aliasing
        spectrum looks ok up to half sample rate, how to recover original smooth signal?
    ENABLE Repeat
        simulating ADC
        would sound tinny, like 8bit video game
        spectrum goes farther out (always to half sample rate, which is now 10x before)
        the spike we want plus higher-frequency crap
    ENABLE Low Pass
        works even at 400kHz when things look like crap
        linear: all of these operations are linear
        Fourier says that any signal is a sum of sines and cosines of the right amplitude
        they would all pass through this process, up to half the sampling freq



2. sin cosine math
    Show 2 cosines
        Change #2
        Change phase
    ENABLE Add
        time plots looks complicated,
        but frequency plot is just the sum (Fourier Transforms are linear)
        change phase; Fourier plots don't show phase
    ENABLE Multiply
        time plot looks equally complicated,
        but frequency plot is now the sum and difference FREQUENCIES
            Board:
                cos(a±b) = cos(a)cos(b) ∓ sin(a)sin(b)  (opposite sign)
                2 cos(a)cos(b) = cos(a+b) + cos(a-b)
                cos(a)cos(b) = 1/2 cos(a+b) + 1/2  cos(a-b)
                similar for sin*sin and sin*cos: always sum and difference freq
                Probably don't waste time writing:
                2 cos(a)cos(b) = cos(a+b) + cos(a-b)
                2 sin(a)sin(b) = cos(a+b) - cos(a-b)
                2 sin(a)cos(b) = sin(a+b) + sin(a-b)  (maybe this one becuase we multiply by cos)
                2 cos(a)sin(b) = sin(a+b) - sin(a-b)
        Try 200 and 10
        sum looks like the sum, mult looks like the product


3. complex exponential math
    Board:
        e^(iθ)  = cos(θ) + i sin(θ)
            say words about Taylor; even odd, i^2==-1
        e^(-iθ) = cos(θ) - i sin(θ)
            add and subtract equations
        e^(+iθ) + e^(-iθ) =  2 cos(θ)
        e^(+iθ) - e^(-iθ) = 2i sin(θ)
        
        cos(θ) =  1/2 e^(+iθ)  +  1/2 e^(-iθ)
        sin(θ) = -i/2 e^(+iθ)  +  i/2 e^(-iθ)
        
            cosine has equal amounts + and - freq with real coefs
            sine had equal and opposite amounts with imaginary coefs
            a sum of sines and cosines would have a complex amount of e^(+iθ) 
                and the complex conjugage amount of e^(-iθ) (only need to know one complex coef)
    Show 2 complex.
        Re and Im parts
        positive vs negative frequencies in leading vs lagging
        ω = 2πf  T=1/f  [keep]  [i=j discuss]
        cos(ωt) = cos(2πft)  goes through a complete cycle as t -> t + T
        complex Fourier xform is just the (magnitude of the complex) coefficient in front of e^(iωt)
        This is why frequency plots for complex signals go + and - freq
        If your signal is real, the neg coefficients are c.c. of the pos ones.
            We only show magnitude, so it looks sym.
        If your signal is complex, it doesn't need to have equal amounts + and -. e^(iωt) itself is example
        the single spike is SIMPLER than cos/sin, which would have two equal and opposite spikes on this this plot
    ENABLE Add
        Again, linear so Fourier transforms add
    ENABLE Multiply
        simple: just one spike a the sum (no 1/2 sum, 1/2 difference)
        because e^(iαt)e^(iβt) = e^[i(α+β)t]



4. real modulation
    Show white noise
        draw independent samples (from Gaussian) at each sample
        frequency spectrum flat
        real, so spectrum symmetric
    ENABLE Low Pass and modulation cosine
        0 freq modulation to hide it at first
        play with low-pass cutoff and transition width (2 parameters)
        samples are no long independent
        this will be our fake smoothed-digital data or fake audio signal
    ENABLE 1st Multiply
        analog multiplier: like amplifier with electrically-controlled gain (Fast: 915MHz, 5GHz)
        each sin and cosine component of the signal gets multiplied by the fast cos
        each turns into the sum and difference frequency on either side of carrier
        mirrors input spectrum around carrier
        200kHz is good
        there are 4 copies of the original single-sided real spectrum:
            x2 because it's a real signal and negative freqs are always c.c. of pos
            x2 because of the mirroring
        Board:
            want to recover original signal, so will multiply again
            cos(θ)^2 = 1/2 cos(θ+θ) + 1/2 cos(θ-θ) = 1/2 (1+cos(2θ))
            each component of original signal gets multiplied by this, 
            so there is an unchanged + mirrors around 2f
    ENABLE 2nd Multiply
        recover original signal by multiplying again
        have extra stuff
        200kHz is good
        full spectrum was useful because each got shifted both up and down,
            with the stuff around 0 getting added together
        we'll LPF to recover
        see in time domain: bunch of noisy all-positive and all-negative.
        hard to recover if demodulating oscillator is even a little off, even just in phase.
    ENABLE Low Pass
        this is also analog
        signals match
        delay is because processing (esp filtering) introduces some digital lag because it needs to accumulate
        went back to one-sided Fourier plot for comparison
        200kHz modulation freq to best illustrate
        50kHz modulation with a narrower LPF cutoff to adjust phase. 0 at 5-sample delay, -1 at 10-sample delay 


5. complex modulation
    Show white noise
        same, but independent draws for Re and Im
        non symmetric spectrum
    ENABLE Low Pass and modulation exp
        0 freq  modulation to hide it
        Re and Im are independent
        spectrum is not symmetric around zero
    ENABLE 1st Multiply
        each component just gets *e^(iωt) which shifts it up
        e^(iαt)e^(iωt) = e^[i(α+ω)t]
        we'll demodulate just by *e^(-iωt)   (the complex conjugate)
    ENABLE 2st Multiply
        that's it
    ENABLE Low Pass
        not needed, but keeps delay the same
        signals match
        50kHz modulation with a narrower LPF cutoff to adjust phase. 0 at 5-sample delay, -1 at 10-sample delay 



6. zero IF modulation
    "Best or worst of both worlds."
    Complex signal -> real RF -> complex signal : just shifting non-symmetric spectrum around.
    Show white noise
    ENABLE Low Pass right away
        so far everything looks the same
    ENABLE 1st Multiplies and Subtract
        what's going on here? Similar to complex mod *e^(iωt), but only keep real part
        z = x+iy
        z*e^(iωt) = (x+iy)*(cos + i sin) = x*cos - y*sin + i(stuff we ignore)
        s = x*cos - y*sin   [just the Re part.]
         Called "pass-band signal": It exists in a band around the RF center freq of same bandwidth
        x(t) is called I(t) for "in phase" component
        y(t) is called Q(t) for the "quadrature" component; name for quarter-cycle, or 90° out of phase
        the real numbers and multiplies are all analog (Fast: 915MHz, 5GHz)
        spectrum is original spectrum shifted up around modulation frequency... 
        ...plus it's mirror because RF is real
        how to recover?  Similar to complex demod *e^(-iωt)
    ENABLE 2nd Multiplies and GUI
        analog demodulating multiplies: take RF s and multiply by cos or sin
        with trig identities....
        s * 2cos(ωt)  = x + x cos(2ωt) - y cos(2ωt)
        s * 2sin(-ωt) = y - y cos(2ωt) - x sin(2ωt)
        see 0 freq along with 2ω terms.
        (extra factor of 2 is just to recover original amplitude; power spread out; we filtered it away
        0 freq are each symmetric, but we're going to form z = x+iy next
    ENABLE Low Pass and Float To Complex
        analog low-pass filter in real RF systems
        analog-to-digital convert both low-passed signals
        SDR sends both samples over USB, which then get interpreted as a single complex number z = x+iy
        matches
        spectrum is not symmetric (Fourier of complex) but matches original.
        50kHz modulation with a narrower LPF cutoff to adjust phase. 0 at 5-sample delay, -1 at 10-sample delay 
    Discussion:
        how much data do we get from SDR?
        If sample rate is 1MHz, we get 2 numbers every MHz. (1 complex Msps)
        The bandwidth is just 1MHz. For complex stuff, 
            we don't have to worry about Nyquist's factor of 2. Factor of 2 comes in Re & Im samples.
        If we center on 915MHz, we get a shifted copy of the spectrum 
            from 914.5 to 915.5 (1MHz wide) no redundant mirrored stuff.
        But we sometimes have to be careful about what this means and where it came from.
        If my demod frequency is a bit off from my modulation frequency, it'll be as if I now have
            not z = x+iy   but   z * e^[i(ω1-w2)t)  
            and as long as ω1-w2 is small (within band), can correct it in software
        Hardware issues:
            cos and sin must be exactly 90 degrees out of phase no matter the freq
               and gains must be precisely the same
               and local oscillator itself leaks into signal, sometimes via antenna or input amp
               and low frequencies are generally noisier (power line hum 1/f noise)
            this leads to images and things like the 0Hz spike in actual SDRs, 
            which is why we tune off-center or put 0Hz between channels if there is room.
            We can later shift down with complex exp.
         








Euler and basic trig identities:
    Sum and difference:
    sin(a±b) = sin(a)cos(b) ± cos(a)sin(b)
    cos(a±b) = cos(a)cos(b) ∓ sin(a)sin(b)  (opposite sign)
    
    Product:
    2 cos(a)cos(b) = cos(a+b) + cos(a-b)
    2 sin(a)sin(b) = cos(a+b) - cos(a-b)
    2 sin(a)cos(b) = sin(a+b) + sin(a-b)
    2 cos(a)sin(b) = sin(a+b) - sin(a-b)
    Lesson: multiplying two sinusoids gives sinusoids at the sum and difference frequencies
    
    Euler's Identities:
    e^(iθ)  = cos(θ) + i sin(θ)
    
    
    e^x =  = 1 + x + x^2/2! + x^3/3! + x^4/4! + x^5/5! + x^6/6! + x^7/7! +  + ...
    sin(x) =     x          - x^3/3!          + x^5/5!          - x^7/7! + ...
    cos(x) = 1     - x^2/2!          + x^4/4!          - x^6/6!          + ...
    
    
    e^(iθ)  = cos(θ) + i sin(θ)
    e^(-iθ) = cos(θ) - i sin(θ)
    e^(+iθ) + e^(-iθ) =  2 cos(θ)
    e^(+iθ) - e^(-iθ) = 2i sin(θ)
    
    e^(iα)e^(iβ) = e^[i(α+β)]
    Lesson: Multiplying complex exponentials together is easier than trig.
    You can prove many of those trig IDs by using Euler's formula on this and equating real and imaginary parts.
    
Fourier Analysis:
    Demo: add sinusoids to get square wave and triangle wave
    Fourier analysis of a real signal is just the coefficients of the sine and cosine components

Complex Fourier Analysis:
    In some sense easier, just pick off the (complex) coefficients of e^(iωt) where ω can have either sign.

Basic sampling and Nyquist Limit:
    You might think 10x, but 2x is good enough
    Demo: 1MHz sampling, change frequency of signal, watch in time and frequency (comment on negative frequencies here?)
    Demo: Orig, 1-in-N down-sample, then up-sample. See aliasing.


Real modulation

Complex modulation

Zero IF modulation


