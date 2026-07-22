### A Pluto.jl notebook ###
# v0.20.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ e7954fc5-457a-4cce-bcbc-646145877a2e
using SegyIO

# ╔═╡ fc9690f6-1f7a-4d00-99c0-3586b3cff97a
using Plots

# ╔═╡ cc729221-ee7c-4182-9aa4-2c627f42b6c6
using PlutoUI

# ╔═╡ 75eb87b0-64f0-4dd0-b8cc-ee433b96b58e
# FFTW: Fastest Fourier Transform in the West
using FFTW

# ╔═╡ 9aabd1f3-ee90-4748-a5f6-c4c67927ee6d
# neaded for "mean" function
using Statistics

# ╔═╡ 5f41d447-299f-4c56-ac0b-f781be29d646
using SpecialFunctions

# ╔═╡ 2fbdf9ba-e982-40d4-838f-53d5b498c017
md"""
# First look at seismic reflection data

In this **Pluto** notebook, we will start analyzing seismic reflection data as geophysicists might.
"""

# ╔═╡ 34b6289c-4b65-4519-8496-0377941f64aa
md"""
## Downloading data 
"""

# ╔═╡ b015f36d-8859-4ab4-a26e-09a8b8ba5a97
if !isfile("class.sgy")
    download("https://www.dropbox.com/scl/fi/gbdanq178owytm7m3h3tu/class2021.sgy?rlkey=510ivt8ywqul0hx85h8tserma&dl=0", "class.sgy")
end

# ╔═╡ 7c7d3ba2-113a-448e-9e80-7000602cb649
md"""
!!! warning

    The file `class.sgy` is large (5.4 Gb) and make take some time to download.
"""

# ╔═╡ 55f0083b-c457-4ee0-b826-fa994f5611a1
md"""
## Reading data in SEGY format

Seismic reflection data are commonly stored in [SEGY format](http://seg.org/publications/tech-stand/seg_y_rev1.pdf). 

The SEGY format is trace-oriented. A SEGY file contains supplemental information in trace headers. 
"""

# ╔═╡ 0f9b62e2-b79d-4043-9d5a-8495df1e360f
seismic = segy_read("class.sgy");

# ╔═╡ 769e3079-b4ab-4450-b13c-07d6ae9d29f7
typeof(seismic)

# ╔═╡ d731add8-fbb7-46c0-815a-961aa257d746
(nt, ntraces) = size(seismic.data)

# ╔═╡ 38c7a939-2dc7-4a73-acb8-7c3f6c01b0fd
md"""
The data portion of the dataset contains 697,761 traces, each containing 2001 samples.
"""

# ╔═╡ 459eeaa8-536e-4136-9ad9-c7f5ebf998d4
typeof(seismic.data)

# ╔═╡ 9a981b8d-ab3a-4498-b798-9c2ecbe2b03c
md"""
The dataset is a matrix of type `IBMFloat32` (32-bit floating-point numbers in IBM format).

However, this matrix represents a three-dimensional cube. We must invoke inline and crossline counters stored in trace headers to extract the cube.
"""

# ╔═╡ 7b260c3a-dcce-43b1-9550-41c9a52802a0
# inline counter
iline = get_header(seismic,:Inline3D)

# ╔═╡ 9e38f746-87b4-4fa2-b13b-fa7198ffc7f5
# crossline counter
xline = get_header(seismic,:Crossline3D)

# ╔═╡ 56d71912-2cc5-4a67-a83c-6b238f7602d3
ilines = maximum(iline) - minimum(iline) + 1

# ╔═╡ 4122c4fd-3ab8-489b-9ee1-f52499bd9c78
xlines = maximum(xline) - minimum(xline) + 1

# ╔═╡ af17853f-84e5-4d11-a1b9-9fadfb9ae178
xlines * ilines == ntraces

# ╔═╡ 00defa7f-81f1-4009-9067-2b7f53fbcfd1
# create a three-dimensional cube and convert data to standard Float32
cube = reshape(Float32.(seismic.data)/1000, (nt, xlines, ilines));

# ╔═╡ b29f084d-d7b6-4d28-b9cb-a5446904fba6
md"""
## Displaying seismic data
"""

# ╔═╡ a4dda64b-4693-4ef7-9b12-840d8575435b
# time sampling in s
dt = get_header(seismic, :dt)[1] * 1e-6

# ╔═╡ 8de22f97-6aa9-4a4e-b892-d935b4be0bb1
time = range(start=0, length=nt, step=dt)

# ╔═╡ d717c0cf-08a4-4a6d-8d35-21fcb2549c04
il = range(start=minimum(iline), stop=maximum(iline))

# ╔═╡ 7bdf1d39-e10e-4dd3-b281-79e33794732b
xl = range(start=minimum(xline), stop=maximum(xline))

# ╔═╡ edcd372d-383b-4259-bb30-840daa7bad81
plot_xline(line; maxnt=nt) = heatmap(xl, time[1:maxnt], cube[1:maxnt,:,line], 
	yflip=:true, cmap=:grays, clim=(-5, 5),
	xlabel="Crossline", ylabel="Time (s)", title="Inline $(line+il[1]-1)")

# ╔═╡ 274e2a4e-823c-40a5-a60f-b1303325bbff
@bind slice Slider(1:ilines, show_value=true)

# ╔═╡ 85fbf6de-995a-4261-8656-3651e7336b02
md"""
## Using color

It is appropriate to use grayscale for the default display of seismic data. However, using color can help highlight features of interest and assist seismic interpretation. A good color table should respect the character of seismic signals.
"""

# ╔═╡ 89eaa52c-0102-4b73-b5e4-f4de2f778676
plot_xline(line, colormap=:grays; maxnt=nt) = heatmap(xl, time[1:maxnt],
	cube[1:maxnt,:,line], yflip=:true, cmap=colormap, clim=(-5, 5),
	xlabel="Crossline", ylabel="Time (s)", 
	title="Inline $(line+il[1]-1), Colormap=$(colormap)")

# ╔═╡ be1352ae-10cd-4323-9c9b-0866d336aad3
plot_xline(200)

# ╔═╡ b7e3d648-7866-450c-b7aa-fa834e29414e
plot_xline(slice, maxnt=1001)

# ╔═╡ dbcc9f51-db4c-4c25-aef7-2548f4e57d96
plot_xline(slice, :viridis)

# ╔═╡ abf38f8f-889e-49d6-ac2f-6b6729c2847d
plot_xline(slice, :balance)

# ╔═╡ ebf3b5d4-55f8-4820-bf18-f6e89a0088cd
md"""
!!! important
    ### Task 1

    Examine different color schemes at **[https://docs.juliaplots.org/latest/generated/colorschemes/](https://docs.juliaplots.org/latest/generated/colorschemes/)**, try them for displaying the seismic section above, and choose your favorite scheme. 

    **Extra Credit** for creating your own color scheme.
"""

# ╔═╡ 0304a92f-c33f-42fe-92d2-8a6c62020d7d
md"""
### Justification for the chosen color scheme

*add your answer here*
"""

# ╔═╡ ccc8afe4-2a9b-4d8b-905a-183b04afe38b
md"""
Naturally, we can slice the seismic cube not only in crosslines but also in other directions, such as time slices.
"""

# ╔═╡ c5a0214d-13b3-46da-85ec-f9375be5dcd8
plot_timeslice(it, colormap=:grays) = heatmap(il, xl, cube[it,:,:], 
	cmap=colormap, clim=(-5, 5),
	xlabel="Crossline", ylabel="Inline", 
	title="Time Slice at $(time[it]) s, Colormap=$(colormap)")

# ╔═╡ 0e3aa3c5-8e17-4f06-adae-d83e62a5232a
plot_timeslice(450)

# ╔═╡ 574df3af-3a35-4065-b759-dd223ab0c7d7
md"""
## Displaying one trace

To better understand the character of seismic data, let us extract one trace from the middle of the cube and examine it closer.
"""

# ╔═╡ b2388809-f2d8-480c-9429-06874b7a12d0
trace = cube[:,250,500];

# ╔═╡ 80643060-df5f-4aae-9766-d1471d07db03
plot(time, trace, title="Seismic Trace", xlabel="Time (s)", ylabel="Amplitude")

# ╔═╡ 43319893-b38f-4e10-849b-f4082b7ddac6
plot_trace(mint,maxt) = plot(time[mint:maxt], trace[mint:maxt], label=:none,
	fillrange=(zeros(maxt-mint+1),max.(zeros(maxt-maxt+1),trace[mint:maxt])),
	title="Seismic Trace", xlabel="Time (s)", ylabel="Amplitude")

# ╔═╡ 0917bb03-1731-4717-b69c-8a4a10d992cc
plot_trace(351,751)

# ╔═╡ 0ade129c-2bdd-4fe4-a732-fa4533fa4573
plot_trace(1351,1751)

# ╔═╡ ff34c68d-8d7a-469a-8737-fa6b9c416f3f
md"""
What can we say about the character of the signal along a seismic trace?

First, it is oscillatory: the negative and positive parts follow each other. Second, the apparent period of oscillations varies with time. 
"""

# ╔═╡ be5e4114-b237-48c2-85ba-71ae11f93913
md"""
## Frequency domain

To fully understand the frequency content of seismic data, we can apply the Fourier transform to move the data from the time domain to the frequency domain.
"""

# ╔═╡ a682adf8-e319-401d-ae1c-db75e7665496
fourier = rfft(trace);

# ╔═╡ ddcbc280-9d73-4be8-b05e-6de16226b761
# frequency range
freq = rfftfreq(nt, 1/dt);

# ╔═╡ 1c4115ef-2274-44fb-a4ce-62c4932a13e6
plot(freq, real(fourier), title="Real Part of the Fourier Transform", 
	 label=:none, xlabel="Frequency (Hz)")

# ╔═╡ 94f0cefa-2c5d-413a-8ba0-22954ed9b246
plot(freq, imag(fourier), title="Imaginary Part of the Fourier Transform", 
	 label=:none, xlabel="Frequency (Hz)")

# ╔═╡ c0ece7eb-ea7d-406e-9bfe-4861c03ddc27
md"""
Theoretically, the Fourier transform is defined in the continuous world as a transformation from a time-domain function $f(t)$ to its frequency-domain counterpart $F(\omega)$ as follows:

$$F(\omega) = \int\limits_{-\infty}^{\infty} f(t)\,e^{-i\omega\,t}\,dt\;.$$

In the discrete world of digital signals, the corresponding transformation is given by DFT (Discrete Fourier Transform):

$$F_k = \sum\limits_{n=0}^{N-1} f_n\,e^{-i\omega_k\,n\,\Delta t}\;,$$

DFT transforms a vector of values $f_n$, representing a function regularly sampled in time $f_n = f(n\,\Delta t)$ to a vector $F_n$ representing the Fourier transform regularly sampled at

$$\omega_k = \frac{2\pi\,k}{N\,\Delta t} \;,\quad \mbox{for} \, k=0,1,2,\cdots,N-1\;.$$
"""

# ╔═╡ ffb2e17b-93d7-4d8a-9283-f91bde20d0e3
md"""
The frequency $f_k = \displaystyle \frac{\omega_k}{2\pi}$ is measured in hertz (cycles per second). The maximum frequency is $1/\Delta t$.
"""

# ╔═╡ 2e150269-0fdd-49e1-8d67-7ca770e29ff6
plot(freq, abs.(fourier), title="Fourier Spectrum", linewidth=2, color=:red,
	 label=:none, xlabel="Frequency (Hz)")

# ╔═╡ 0c70270e-f539-4abe-aebc-1bcb9bb86452
md"""
The Fourier transform produces a complex-valued signal. 

The spectrum is its absolute value $S(\omega) = |F(\omega)|$.

Let us compute the average spectrum of the whole dataset.
"""

# ╔═╡ 4f7c49f0-0c8e-4e46-b925-0cd8a808bc63
# Design a 1-D Fourier transform
ft = plan_rfft(trace);

# ╔═╡ c95d1330-16da-4949-8518-e1f870ec00c5
# Apply it to every trace in the cube, average across traces
spectrum = mean(abs.(mapslices(t -> ft * t, cube, dims=1)), dims=(2, 3));

# ╔═╡ 50c4d2e8-b93e-4e5e-be91-543e46dbec11
# peak frequency
peak = freq[argmax(spectrum)]

# ╔═╡ 8c29863d-75c5-4460-b6bb-c9dedf63ffb9
# centroid frequency
cent = sum(spectrum .* freq)/sum(spectrum)

# ╔═╡ 85d9d2f4-b969-46cb-b6fb-e090c7ce64ed
begin
	plot(freq, dropdims(spectrum, dims=(2, 3)), title="Average Spectrum", 
		 linewidth=2, label=:none, xlabel="Frequency (Hz)")
	plot!([peak, peak], [0, 52], linestyle=:dash, label="peak frequency")
	plot!([cent, cent], [0, 52], linestyle=:dash, label="centroid frequency")
end

# ╔═╡ d3851a5d-ec2b-49dd-9ca5-2878ec1448ab
md"""
What can we notice from observing the average spectrum?

The data are *band-limited*, containing multiple frequencies inside a frequency band. Frequencies outside of the band are missing. That includes very low frequencies (below about 5 Hz) and high frequencies (above about 90 Hz). The absence of these frequencies limits the ability of seismic reflection data to resolve small subsurface features. 
"""

# ╔═╡ 3db11431-1932-4b6d-96bc-fa425bb9028c
md"""
!!! important

    ### Task 2

    Divide the data into the top 4 seconds and bottom 4 seconds, then compute and plot their corresponding average spectra. What do you notice?
"""

# ╔═╡ e53c0059-e4eb-4e61-8b85-405c8d1ff23e
# window the top 1001 time samples
top = cube[1:1001,:,:];

# ╔═╡ 2ed8e26e-0a99-49d7-86fc-f7d6f7a9e00b
# window the bottom 1001 time samples
# bottom = cube[XXXX];
# uncomment and fix the code above

# ╔═╡ 30065ca1-e437-4208-bc45-c624088c09c5
# Add the average spectrum calculation for top and bot

# ╔═╡ c0ee5397-b70a-46e5-876b-fe30cc26c8ff


# ╔═╡ 553f90fe-2fcf-4299-9a37-154cf0151ea2
md"""
### Summary of comparing the two spectra

*add your answer here*
"""

# ╔═╡ 63c38ba6-4a00-410f-9851-b00c5cbe467c
md"""
## 2-D Fourier transform

The concept of transforming data extends to multiple dimensions. 

The 2-D Fourier transform takes the data from time and space $u(t,x)$ to frequency and wavenumber (spatial frequency) $U(\omega,k)$. In the continuous world, the transformation is defined by the following 2-D integral:

$$U(\omega,k) = \iint\limits_{-\infty}^{\infty} u(t,x)\,e^{-i\omega\,t+ik\,x}\,dt\,dx\;.$$
"""

# ╔═╡ 2da5d684-92c2-4686-ad12-e1ac6029e245
section = cube[:,:,slice];

# ╔═╡ ff71f09d-f2fa-4a09-9383-a9d05555a251
# Apply 2-D real-to-complex Fourier transform
fourier2 = mapslices(x -> fftshift(fft(x)), 
	       mapslices(t -> ft * t, section, dims=1), dims=2);

# ╔═╡ 69217abd-5d47-4ec6-9681-e2e832ae7abb
# the trace spacing is 25 meters
dx = 0.025

# ╔═╡ 39f60659-52f4-4720-8208-b2d36122b1f6
# wavenumbers
waven = fftshift(fftfreq(xlines, 1/dx));

# ╔═╡ 8f500daa-d2d6-45e1-8642-af6f146eaa18
begin
	heatmap(waven, freq, abs.(fourier2), clim=(0,5000),
	        title="2-D Fourier Spectrum", yflip=:true,
		    xlabel="Wavenumber (1/km)", ylabel="Frequency (Hz)")
	plot!([-20, 0, 20], [40, 0, 40], color=:white, linestyle=:dash, label=:none)
	plot!([-20, 20], [75, 75], color=:white, linestyle=:dash, label=:none)
end

# ╔═╡ 627b2c87-2e8d-42e3-8c73-7ca972506927
md"""
If each point in a 1-D Fourier transform corresponds to a particular frequency, each point in a 2-D Fourier transform corresponds to a plane wave $w(t-x/v)$ with the apparent velocity $v=\omega/k$.

Nearly horizontal events in the original domain ($v=\infty$) appear near vertical in the 2-D Fourier domain ($k=0$), and vice versa. 

Looking at the 2-D spectrum of a seismic section, we notice a band of slopes in addition to the previously established band of temporal frequencies. The slope limitation indicates a difficulty for seismic reflection data recorded at the Earth's surface in resolving events with near-vertical slopes.
"""

# ╔═╡ 1b31ec1d-8bfb-4427-82af-d7a4532e054c
md"""
## "Seismic Lena"

To better understand the limitations of seismic data, we will try to reproduce an example [originally suggested](https://library.seg.org/doi/10.1190/1.1438642) by Chris Liner, a famous geophysicist. 

![](https://wehco.media.clients.ellingtoncms.com/img/photos/2015/01/23/resized_99261-nwpliner0201-color-cover_15-19233_t800.JPG?90232451fbcadccc64a17de7521d859a8f88077d)

See also the follow-up paper by David Monk.

* Monk, D., 2002. [Lena: A seismic model](https://library.seg.org/doi/10.1190/1.1481249). The Leading Edge, 21(5), pp.438-444.
"""

# ╔═╡ 21bd1cbe-0e5c-42ce-8d48-143d5bb01961
download("https://fomel.com/data/imgs/lena.img", "lena.img")

# ╔═╡ b1197e65-d83e-4f08-abd6-1ade879e1c3f
# binary data
lena = Array{UInt8}(undef, 512, 513);

# ╔═╡ 96508850-dbe8-49ec-bd32-dda1e2e09454
read!("lena.img", lena);

# ╔═╡ cd8fab47-5be8-4537-bcb5-4762093a344e
heatmap(lena', yflip=:true, color=:grays, title="Lena",
	    showaxis=:false, aspect_ratio=1, grid=:false, colorbar=:false)

# ╔═╡ f1ed700f-b570-4728-8c25-35b3312dc898
md"""
"Lena" is an image that played an important historical role in the development of image analysis algorithms. Its usage is no longer recommended.

[https://en.wikipedia.org/wiki/Lenna](https://en.wikipedia.org/wiki/Lenna)
"""

# ╔═╡ f44963fe-9104-4c7c-8cb6-88b97a2921e5
# transpose and convert to floating point
flena = Float32.(lena');

# ╔═╡ 7bf14232-dd66-4b09-841b-6c92f63b79f9
lena_f = mapslices(x -> fftshift(fft(x)), 
	               mapslices(t -> rfft(t), 
						     flena, dims=1), dims=2);

# ╔═╡ fbba9c51-9988-41c9-92f7-192eee76ec8f
lena_freq = rfftfreq(513, 1/dt);

# ╔═╡ 19d146a6-6806-4b85-9f81-526176139b81
lena_waven = fftshift(fftfreq(512, 1/dx));

# ╔═╡ a378ffd3-f537-4476-8da4-534883592138
heatmap(lena_waven, lena_freq, abs.(lena_f)/1f5, 
	    yflip=:true, title="2-D Lena Spectrum", clim=(0, 1),
	    xlabel="horizontal frequency", ylabel="vertical frequency")

# ╔═╡ c076af4e-e20f-43fd-a6a2-b520bfb9b2f1
# bandpass filtering by windowing in the Fourier domain
function bandpass(ft, fs; f1=:none, f2=:none, a=1)
    nf, n2 = size(ft)
	F = copy(ft)
	for i2 in 1:n2
	    for i in 1:nf
    	    f = fs[i]
			if f1 != :none
            	F[i,i2] *= 1 - (erf(a*(f1 + f)) + erf(a*(f1 - f)))/2
			end
			if f2 != :none
			 	F[i,i2] *= (erf(a*(f2 + f)) + erf(a*(f2 - f)))/2
			end
        end
    end
    return F
end

# ╔═╡ 463c2eb7-6edc-407a-8c72-207b1f3837ed
lena_filt = bandpass(lena_f, lena_freq, f1=10, f2=40);

# ╔═╡ 392c0434-f66e-48b4-b603-55012e22bcca
heatmap(lena_waven, lena_freq, abs.(lena_filt)/1f5, 
	    yflip=:true, title="2-D Lena Spectrum (Band-limited)", clim=(0, 1),
	    xlabel="horizontal frequency", ylabel="vertical frequency")

# ╔═╡ ba12c76e-777f-427b-8ae2-b9e0a9db0bd3
# slope filtering by windowing in the Fourier domain
function dipfilter(ft, fs, ks, vel; a=1)
    nf, n2 = size(ft)
	F = copy(ft)
	for i2 in 1:n2
		f1 = abs(ks[i2])*vel
	    for i in 1:nf
    	    f = fs[i]
			F[i,i2] *= 1 - (erf(a*(f1 + f)) + erf(a*(f1 - f)))/2
        end
    end
    return F
end

# ╔═╡ fed79755-34d0-47d6-877e-4d860b728831
lena_filt2 = dipfilter(lena_filt, lena_freq, lena_waven, 1);

# ╔═╡ 706557c9-f24d-409a-a690-f5d96b94621b
heatmap(lena_waven, lena_freq, abs.(lena_filt2)/1f5, 
	    yflip=:true, title="2-D Lena Spectrum (Bandlimited and Diplimited)", 
	    clim=(0, 1), xlabel="horizontal frequency", ylabel="vertical frequency")

# ╔═╡ ec6bf25d-d68e-4016-b979-036d19e08df9
# Apply the inverse 2-D Fourier transform
seismic_lena = mapslices(t -> irfft(t, 513), 
	               mapslices(x -> ifft(ifftshift(x)), 
						     lena_filt2, dims=2), dims=1);

# ╔═╡ 3bdd0f5c-da7c-41bf-8850-04995cc23561
heatmap(seismic_lena, yflip=:true, color=:grays, title="Seismic Lena",
	    showaxis=:false, aspect_ratio=1, grid=:false, colorbar=:false)

# ╔═╡ 26231eb0-32eb-471e-b707-4f3f8282e3de
md"""
We can see that band-limited and dip-limited data retain some original information but lose resolution and absolute scale.

Because of the dip limitation, the information content is limited to near horizontal edges while missing nearly vertical edges.
"""

# ╔═╡ 2bd6c464-2d08-4452-bee1-6e37d4a27917
md"""
## Fundamentals of seismic wave propagation

To understand the theory of seismic wave propagation, let us first consider the propagation of waves in fluids. To derive the corresponding mathematical equations, we can start with Newton's second law of motion: mass times acceleration is equal to the force causing the motion:

$$\displaystyle m\,\frac{d^2\mathbf{x}}{dt^2} = \mathbf{F}\;.$$

Inside a fluid, the mass of a small volume with dimensions $\Delta x$, $\Delta y$, and $\Delta z$ is the density $\rho$ times the volume. Suppose for simplicity that the displacement happens only in the vertical direction $z$ and is caused by the change in pressure $P$ between the top and the bottom of the volume.

$$\displaystyle \rho\,\Delta x\,\Delta y\,\Delta z\,\frac{d^2 u_z}{dt^2} = \left[P(z)-P(z+\Delta z)\right]\,\Delta x\,\Delta y\;.$$

Dividing both sides of the equation by the unit volume and setting $\Delta z$ to the infinitely small limit, we arrive at the following differential equation, which expresses Newton's second law for the case of fluid motion in a continuous media:

$$\displaystyle  \rho\,\frac{d^2 u_z}{dt^2} = - \frac{d P}{d z}.$$
"""

# ╔═╡ b20e5678-d42a-439b-94ca-3a11550482d6
md"""
While changes in pressure cause fluid motion, the motion causes changes in pressure. Unlike Newton's law, our second equation is not a fundamental law of physics but an approximation. The amount of fluid leaving a unit volume is $\left[u_z(z+\Delta z)-u_z(z)\right]\,\Delta x\,\Delta y$, or, per unit volume and in the limit of infinitely small increment, the derivative $d u_z/d z$. We will assume that the loss of fluid causes a drop in pressure according to

$$P = \displaystyle - K\,\frac{d u_z}{d z}\;,$$

where $K$ is a material property (known as *bulk modulus* and having the physical units of pressure). Putting the two equations together leads to the second differential equation describing the wave motion:

$$\displaystyle  \rho\,\frac{d^2 u_z}{dt^2} = \frac{d}{d z} \left(K\,\frac{d u_z}{d z}\right)\;.$$
"""

# ╔═╡ eeeddad1-4ce6-40a6-a336-abe84a4a8576
md"""
In a homogeneous medium, a plane wave propagating with a constant velocity

$$u_z(z,t) = \displaystyle f\left(t-\frac{z}{v}\right)$$

will satisfy the wave equation provided that 

$$v^2 = \displaystyle \frac{K}{\rho}\;.$$

The corresponding pressure is

$$P(z,t) = \displaystyle \frac{K}{v} f'\left(t-\frac{z}{v}\right) = v\,\rho\,\frac{d u_z}{d t}\;.$$

The product of velocity and density $v\,\rho$, which appears in the equation above, is known as *acoustic impedance*: 

$$I = v\,\rho\;.$$

The acoustic impedance is the proportionality coefficient between pressure and the velocity of fluid motion for the vertically propagating plane wave.
"""

# ╔═╡ 85e5ec2a-c062-42ff-9b11-549a8f971988
md"""
## Fundamentals of seismic reflection

Let us now consider a horizontal interface between two layers with different physical properties. As an incident wave hits the interface, it splits into the reflected and transmitted parts. The reflected wave travels back into the first layer, and the transmitted part enters the second layer. Across the interface, the motion velocity $v = d u_z/dt$ needs to remain continuous:

$$v_I - v_R = v_T\;,$$

where $v_I$, $v_R$, and $v_T$ represent the incident, reflected, and transmitted wavefields. The pressure across the interface is also continuous:

$$P_I + P_R =  P_T$$

or, equivalently,

$$I_1\,v_I + I_1\,v_R =  I_2\,v_T\;,$$

where $I_1$ and $I_2$ represent acoustic impedances in the two layers. Putting the equations for velocity and pressure continuity together, we can derive an expression for the reflection coefficient

$$r = \displaystyle \frac{v_R}{v_I} = \frac{I_2 - I_1}{I_2 + I_1}\;.$$
"""

# ╔═╡ af1a920e-14a2-4011-b1bd-881aefeaab65
md"""
Note the many strong assumptions that we made in this derivation:

1. Acoustic media.
2. Plane waves.
3. A plane interface.
4. Normal incidence.

Nevertheless, the reflection coefficient formula is helpful because it indicates the properties of the media to which seismic data are sensitive.
"""

# ╔═╡ 2ed499a1-9d7c-43f4-b6ab-47a615fdfd12
md"""
!!! note
    The seismic reflection method aims to image the subsurface reflectivity, which indicates contrasts in acoustic impedance at geological interfaces.
"""

# ╔═╡ d236ccc9-70a0-4e1f-9de9-967c36cec91c
md"""
## Marmousi-2

Marmousi is a famous synthetic Earth model inspired by geological formations offshore West Africa. The geophysical community used it in numerous computational experiments.

The original Marmousi model from the 1990s contained only seismic velocity. In the 2000s, it was extended to a larger grid and additional physical properties, such as density.

* Versteeg, R., 1994. [The Marmousi experience: Velocity model determination on a synthetic complex data set](https://library.seg.org/doi/abs/10.1190/1.1437051). The Leading Edge, 13(9), pp.927-936.
* Martin, G.S., Wiley, R. and Marfurt, K.J., 2006. [Marmousi2: An elastic upgrade for Marmousi](https://library.seg.org/doi/full/10.1190/1.2172306). The Leading Edge, 25(2), pp.156-166.
* [https://wiki.seg.org/wiki/AGL_Elastic_Marmousi](https://wiki.seg.org/wiki/AGL_Elastic_Marmousi)
"""

# ╔═╡ 27ff5aaf-9aa3-4647-89c5-9a7762a762c3
if !isfile("velocity.segy")
    download("https://www.dropbox.com/scl/fi/5uyzl4pabyy0zfpsmihik/MODEL_P-WAVE_VELOCITY_1.25m.segy?rlkey=0dia79gka88cgfuybhnm9t18w&st=u3qn4bml&dl=0", "velocity.segy")
end

# ╔═╡ a7a160f3-3112-48b8-ba50-cbc647d9f41a
if !isfile("density.segy")
    download("https://www.dropbox.com/scl/fi/7c9lom68u5qqyghlxshlr/MODEL_DENSITY_1.25m.segy?rlkey=32on42qzvs4o9sampqigupwhl&st=v46xteq4&dl=0", "density.segy")
end

# ╔═╡ 924b8f3d-67ef-487d-ac5c-e1dc478b40f0
velocity = segy_read("velocity.segy");

# ╔═╡ 73b4ae8e-ea37-4f80-8632-b23d97876905
density = segy_read("density.segy");

# ╔═╡ 95b84958-19fb-4755-83c5-00171aa15287
(nz, nx) = size(velocity.data)

# ╔═╡ 406a1a3e-fd12-4190-8842-b858727b1173
# grid size
dz = 0.0025

# ╔═╡ 067527ac-621f-4f5e-8d42-c3ebf07cde60
depth = range(start=0, length=nz, step=dz)

# ╔═╡ f654d24a-2fff-4a51-b919-2d0abfa0e014
lateral = range(start=0, length=nx, step=dz)

# ╔═╡ 0998f8ef-89e9-4532-9bfa-343820b03605
vel = Float32.(velocity.data)/1000; # velocity in km/s

# ╔═╡ 559de934-506a-4e8d-a97b-d785c3ab1bc1
heatmap(lateral, depth, vel, yflip=:true,
	xlabel="Lateral (km)", ylabel="Depth (km)", cmap=:viridis,
	title="Marmousi-2 Velocity (km/s)")

# ╔═╡ f7d9a74f-7ec8-4f2f-ba65-e1d51a8d2c83
heatmap(lateral, depth, Float32.(density.data), yflip=:true,
	xlabel="Lateral (km)", ylabel="Depth (km)", 
	title="Marmousi-2 Density (g/cm^3)")

# ╔═╡ 83620794-1d4b-4464-a913-4fdefe999a52
impedance = Float32.(density.data) .* vel;

# ╔═╡ 84c19b7b-2dd8-4622-b4ea-1e59f508af08
heatmap(lateral, depth, impedance, yflip=:true,
	xlabel="Lateral (km)", ylabel="Depth (km)", cmap=:viridis,
	title="Marmousi-2 Impedance")

# ╔═╡ 1b320ce6-9ed4-485f-83ca-9128736c62bc
md"""
!!! important

    ## Task3
    Convert the Marmousi-2 acoustic impedance to reflectivity and display the result.

    **Extra credit** for removing low and high frequencies using bandpass filtering.
"""

# ╔═╡ 002e4c8a-bb93-40d3-829c-9ceb297ce7e5
# convert acoustic impedance to reflectivity
function ai2refl(imp)
	nz, nx = size(imp)
	ref = similar(imp)
	ref0 = zero(eltype(imp))
	for ix in 1:nx
		imp1 = imp[1, ix]
		for iz in 1:nz-1
			imp2 = imp[iz+1, ix] 
			ref[iz, ix] = (imp2 - imp1)/(imp2 + imp1)
			imp1 = imp2
		end
		ref[nz, ix] = ref0
	end
	return ref
end		

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SegyIO = "157a0f19-4d44-4de5-a0d0-07e2f0ac4dfa"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
FFTW = "~1.8.1"
Plots = "~1.40.20"
PlutoUI = "~0.7.83"
SegyIO = "~0.8.5"
SpecialFunctions = "~2.4.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "e8bc9ddc6155df8079a49748fe9df352ea66b1a9"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "bbe1079eecf9c9fbb52765193ad2bae27ae09bc8"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.10"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1fa950ebc3e37eccd51c6a8fe1f92f7d86263522"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.7+0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e6c4a6407a949e79a9d3f249bf49e6987c80e01f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.8.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "7a58e45171b63ed4782f2d36fdee8713a469e6e0"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.1.2+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "7de7c78d681078f027389e067864a8d53bd7c3c9"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6866aec60ef98e3164cd8d6855225684207e9dff"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.12+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "9e0fb9e54594c47f278d75063980e43066e26e20"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "f954322d5de03ec630d177cda203dcd92b6be399"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.26"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "6fada551286ab6ea4ca1628cb2de9f166a2ec966"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.26+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "69ffb934a5c5b7e086a0b4fee3427db2556fba6e"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.16+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "c89d196f5ffb64bfbf80985b699ea913b0d2c211"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.6.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1dae3057da6f2b9c857afef03177bbdc7c4afe92"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b7970cef8ae1c990ba0c09cd8bdc1145e006632f"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "22.1.7+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "24390f715ff0795a1c4b912d788f18c52c6abd19"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.11"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58e5ed5e386e156bd93e86b305ebd21ac63d2d04"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "32a4e09c5f29402573d673901778a0e03b0807b9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.6"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "e4a6721aa89e62e5d4217c0b21bd714263779dda"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.46.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "bfe839e9668f0c58367fb62d8757315c0eac8777"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.20"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "144895f6166994730ee7ff8113b981fc360638f1"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+2"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "159d253ab126d5b29230cf53521899bea4ef4648"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SegyIO]]
deps = ["Distributed", "Printf", "Test"]
git-tree-sha1 = "0fc24db28695a80aa59c179372b11165d371188a"
uuid = "157a0f19-4d44-4de5-a0d0-07e2f0ac4dfa"
version = "0.8.5"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "e4d7a1a0edc20af42689ea6f4f3587a2175d50ee"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.12"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "82bee338d650aa515f31866c460cb7e3bcef90b8"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "6258d453843c466d84c17a58732dda5deeb8d3af"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.24.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "af305cc62419f9bd61b6644d19170a4d258c7967"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.7.0"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "58972370b81423fc546c56a60ed1a009450177c3"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.19.0+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "2e59214e017a55cb87474a00fa76035c82ac0e17"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.47.0+2"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "850b06095ee71f0135d644ffd8a52850699581ed"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.3+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "da8c1f6eee04831f14edcfa5dae611d309807e57"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.3.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─2fbdf9ba-e982-40d4-838f-53d5b498c017
# ╟─34b6289c-4b65-4519-8496-0377941f64aa
# ╠═b015f36d-8859-4ab4-a26e-09a8b8ba5a97
# ╟─7c7d3ba2-113a-448e-9e80-7000602cb649
# ╟─55f0083b-c457-4ee0-b826-fa994f5611a1
# ╠═e7954fc5-457a-4cce-bcbc-646145877a2e
# ╠═0f9b62e2-b79d-4043-9d5a-8495df1e360f
# ╠═769e3079-b4ab-4450-b13c-07d6ae9d29f7
# ╠═d731add8-fbb7-46c0-815a-961aa257d746
# ╟─38c7a939-2dc7-4a73-acb8-7c3f6c01b0fd
# ╠═459eeaa8-536e-4136-9ad9-c7f5ebf998d4
# ╟─9a981b8d-ab3a-4498-b798-9c2ecbe2b03c
# ╠═7b260c3a-dcce-43b1-9550-41c9a52802a0
# ╠═9e38f746-87b4-4fa2-b13b-fa7198ffc7f5
# ╠═56d71912-2cc5-4a67-a83c-6b238f7602d3
# ╠═4122c4fd-3ab8-489b-9ee1-f52499bd9c78
# ╠═af17853f-84e5-4d11-a1b9-9fadfb9ae178
# ╠═00defa7f-81f1-4009-9067-2b7f53fbcfd1
# ╟─b29f084d-d7b6-4d28-b9cb-a5446904fba6
# ╠═a4dda64b-4693-4ef7-9b12-840d8575435b
# ╠═8de22f97-6aa9-4a4e-b892-d935b4be0bb1
# ╠═d717c0cf-08a4-4a6d-8d35-21fcb2549c04
# ╠═7bdf1d39-e10e-4dd3-b281-79e33794732b
# ╠═fc9690f6-1f7a-4d00-99c0-3586b3cff97a
# ╠═edcd372d-383b-4259-bb30-840daa7bad81
# ╠═be1352ae-10cd-4323-9c9b-0866d336aad3
# ╠═cc729221-ee7c-4182-9aa4-2c627f42b6c6
# ╠═b7e3d648-7866-450c-b7aa-fa834e29414e
# ╠═274e2a4e-823c-40a5-a60f-b1303325bbff
# ╟─85fbf6de-995a-4261-8656-3651e7336b02
# ╠═89eaa52c-0102-4b73-b5e4-f4de2f778676
# ╠═dbcc9f51-db4c-4c25-aef7-2548f4e57d96
# ╠═abf38f8f-889e-49d6-ac2f-6b6729c2847d
# ╟─ebf3b5d4-55f8-4820-bf18-f6e89a0088cd
# ╠═0304a92f-c33f-42fe-92d2-8a6c62020d7d
# ╟─ccc8afe4-2a9b-4d8b-905a-183b04afe38b
# ╠═c5a0214d-13b3-46da-85ec-f9375be5dcd8
# ╠═0e3aa3c5-8e17-4f06-adae-d83e62a5232a
# ╟─574df3af-3a35-4065-b759-dd223ab0c7d7
# ╠═b2388809-f2d8-480c-9429-06874b7a12d0
# ╠═80643060-df5f-4aae-9766-d1471d07db03
# ╠═43319893-b38f-4e10-849b-f4082b7ddac6
# ╠═0917bb03-1731-4717-b69c-8a4a10d992cc
# ╠═0ade129c-2bdd-4fe4-a732-fa4533fa4573
# ╟─ff34c68d-8d7a-469a-8737-fa6b9c416f3f
# ╟─be5e4114-b237-48c2-85ba-71ae11f93913
# ╠═75eb87b0-64f0-4dd0-b8cc-ee433b96b58e
# ╠═a682adf8-e319-401d-ae1c-db75e7665496
# ╠═ddcbc280-9d73-4be8-b05e-6de16226b761
# ╠═1c4115ef-2274-44fb-a4ce-62c4932a13e6
# ╠═94f0cefa-2c5d-413a-8ba0-22954ed9b246
# ╟─c0ece7eb-ea7d-406e-9bfe-4861c03ddc27
# ╟─ffb2e17b-93d7-4d8a-9283-f91bde20d0e3
# ╠═2e150269-0fdd-49e1-8d67-7ca770e29ff6
# ╟─0c70270e-f539-4abe-aebc-1bcb9bb86452
# ╠═4f7c49f0-0c8e-4e46-b925-0cd8a808bc63
# ╠═9aabd1f3-ee90-4748-a5f6-c4c67927ee6d
# ╠═c95d1330-16da-4949-8518-e1f870ec00c5
# ╠═50c4d2e8-b93e-4e5e-be91-543e46dbec11
# ╠═8c29863d-75c5-4460-b6bb-c9dedf63ffb9
# ╠═85d9d2f4-b969-46cb-b6fb-e090c7ce64ed
# ╟─d3851a5d-ec2b-49dd-9ca5-2878ec1448ab
# ╟─3db11431-1932-4b6d-96bc-fa425bb9028c
# ╠═e53c0059-e4eb-4e61-8b85-405c8d1ff23e
# ╠═2ed8e26e-0a99-49d7-86fc-f7d6f7a9e00b
# ╠═30065ca1-e437-4208-bc45-c624088c09c5
# ╠═c0ee5397-b70a-46e5-876b-fe30cc26c8ff
# ╠═553f90fe-2fcf-4299-9a37-154cf0151ea2
# ╟─63c38ba6-4a00-410f-9851-b00c5cbe467c
# ╠═2da5d684-92c2-4686-ad12-e1ac6029e245
# ╠═ff71f09d-f2fa-4a09-9383-a9d05555a251
# ╠═69217abd-5d47-4ec6-9681-e2e832ae7abb
# ╠═39f60659-52f4-4720-8208-b2d36122b1f6
# ╠═8f500daa-d2d6-45e1-8642-af6f146eaa18
# ╟─627b2c87-2e8d-42e3-8c73-7ca972506927
# ╟─1b31ec1d-8bfb-4427-82af-d7a4532e054c
# ╠═21bd1cbe-0e5c-42ce-8d48-143d5bb01961
# ╠═b1197e65-d83e-4f08-abd6-1ade879e1c3f
# ╠═96508850-dbe8-49ec-bd32-dda1e2e09454
# ╠═cd8fab47-5be8-4537-bcb5-4762093a344e
# ╟─f1ed700f-b570-4728-8c25-35b3312dc898
# ╠═f44963fe-9104-4c7c-8cb6-88b97a2921e5
# ╠═7bf14232-dd66-4b09-841b-6c92f63b79f9
# ╠═fbba9c51-9988-41c9-92f7-192eee76ec8f
# ╠═19d146a6-6806-4b85-9f81-526176139b81
# ╠═a378ffd3-f537-4476-8da4-534883592138
# ╠═5f41d447-299f-4c56-ac0b-f781be29d646
# ╠═c076af4e-e20f-43fd-a6a2-b520bfb9b2f1
# ╠═463c2eb7-6edc-407a-8c72-207b1f3837ed
# ╠═392c0434-f66e-48b4-b603-55012e22bcca
# ╠═ba12c76e-777f-427b-8ae2-b9e0a9db0bd3
# ╠═fed79755-34d0-47d6-877e-4d860b728831
# ╠═706557c9-f24d-409a-a690-f5d96b94621b
# ╠═ec6bf25d-d68e-4016-b979-036d19e08df9
# ╠═3bdd0f5c-da7c-41bf-8850-04995cc23561
# ╟─26231eb0-32eb-471e-b707-4f3f8282e3de
# ╟─2bd6c464-2d08-4452-bee1-6e37d4a27917
# ╟─b20e5678-d42a-439b-94ca-3a11550482d6
# ╟─eeeddad1-4ce6-40a6-a336-abe84a4a8576
# ╟─85e5ec2a-c062-42ff-9b11-549a8f971988
# ╟─af1a920e-14a2-4011-b1bd-881aefeaab65
# ╟─2ed499a1-9d7c-43f4-b6ab-47a615fdfd12
# ╟─d236ccc9-70a0-4e1f-9de9-967c36cec91c
# ╠═27ff5aaf-9aa3-4647-89c5-9a7762a762c3
# ╠═a7a160f3-3112-48b8-ba50-cbc647d9f41a
# ╠═924b8f3d-67ef-487d-ac5c-e1dc478b40f0
# ╠═73b4ae8e-ea37-4f80-8632-b23d97876905
# ╠═95b84958-19fb-4755-83c5-00171aa15287
# ╠═406a1a3e-fd12-4190-8842-b858727b1173
# ╠═067527ac-621f-4f5e-8d42-c3ebf07cde60
# ╠═f654d24a-2fff-4a51-b919-2d0abfa0e014
# ╠═0998f8ef-89e9-4532-9bfa-343820b03605
# ╠═559de934-506a-4e8d-a97b-d785c3ab1bc1
# ╠═f7d9a74f-7ec8-4f2f-ba65-e1d51a8d2c83
# ╠═83620794-1d4b-4464-a913-4fdefe999a52
# ╠═84c19b7b-2dd8-4622-b4ea-1e59f508af08
# ╟─1b320ce6-9ed4-485f-83ca-9128736c62bc
# ╠═002e4c8a-bb93-40d3-829c-9ceb297ce7e5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
