### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> chapter = 4
#> section = 2
#> order = 3
#> image = "https://github.com/quentingageclark-a11y/plutoimages/blob/main/adreg.png?raw=true"
#> title = "First Look"
#> tags = ["module4", "track_julia", "track_material", "Pluto", "PlutoUI"]
#> layout = "layout.jlhtml"


using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
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
download("https://www.dropbox.com/scl/fi/gbdanq178owytm7m3h3tu/class2021.sgy?rlkey=510ivt8ywqul0hx85h8tserma&dl=0", "class.sgy")

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

# ╔═╡ 319b13b9-5536-4ca9-ac3a-8fc19307be7d
download("https://www.dropbox.com/scl/fi/5uyzl4pabyy0zfpsmihik/MODEL_P-WAVE_VELOCITY_1.25m.segy?rlkey=0dia79gka88cgfuybhnm9t18w&st=u3qn4bml&dl=0", "velocity.segy")

# ╔═╡ daf87367-a25a-4ca6-9b94-dfb79eb0c372
download("https://www.dropbox.com/scl/fi/7c9lom68u5qqyghlxshlr/MODEL_DENSITY_1.25m.segy?rlkey=32on42qzvs4o9sampqigupwhl&st=v46xteq4&dl=0", "density.segy")

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
FFTW = "~1.8.0"
Plots = "~1.40.5"
PlutoUI = "~0.7.59"
SegyIO = "~0.8.5"
SpecialFunctions = "~2.4.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "92784a00066c819e68353bcef77ef894184d2cc4"

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
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "f389674c99bfcde17dc57454011aa44d5a260a40"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e16271d212accd09d52ee0ae98956b8a05c4b626"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "17.0.6+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

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

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "082f0c4b70c202c37784ce4bfbc33c9f437685bf"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.5"

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
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

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
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SegyIO]]
deps = ["Distributed", "Printf", "Test"]
git-tree-sha1 = "0fc24db28695a80aa59c179372b11165d371188a"
uuid = "157a0f19-4d44-4de5-a0d0-07e2f0ac4dfa"
version = "0.8.5"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

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

[[deps.TranscodingStreams]]
git-tree-sha1 = "e84b3a11b9bece70d14cce63406bbc79ed3464d2"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.2"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

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
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
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
# ╠═319b13b9-5536-4ca9-ac3a-8fc19307be7d
# ╠═daf87367-a25a-4ca6-9b94-dfb79eb0c372
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
