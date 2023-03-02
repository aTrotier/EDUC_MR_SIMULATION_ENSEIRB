### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ 2e4cbd0e-a92b-11ed-022c-eda67b72483f
begin
	using CairoMakie, PlutoUI, LinearAlgebra,Statistics,Random,NIfTI
end

# ╔═╡ 11bf61dc-aa97-48b6-ab67-9bc0861120bc
md"
# TP IRM SIMULATION
_Aurélien Trotier : Février 2023_

Le signal IRM peut-être représenté grâce à l'équation de Bloch :

$\frac{d\vec{M}}{dt} = \gamma \vec{M} \times \vec{B}_{ext} + \frac{1}{T_1} (M_0 - M_z)\vec{z} - \frac{1}{T_2} \vec{M}_\perp$

Cette équation peut-être résolue dans certaines conditions : 
- Rotation (voir wikipedia)
- Relaxation

Nous utiliserons ici une implémentation en matrice 4x4, ce qui donne pour la relaxation

$$\begin{equation}
A = \begin{bmatrix}
E2 & 0 & 0 & 0\\
0 & E2 & 0 & 0\\
0 & 0 & E1 & 1-E1\\
0 & 0 & 0 & 1
\end{bmatrix}
\end{equation}$$

## 1. Préparation des sous-fonctions

Préparez les fonctions nécessaires à la simulation du signal IRM :

- xrot(angle), yrot(angle) et zrot(angle) : qui retourne une matrice de rotation à partir d'un angle (en radian)
- throt(phi,theta) : si l'angle de rotation est transverse (différent de x/y) (utilisez les fonctions précédentes : $Rth(phi,theta)=Rz(theta)*Rx(phi)*Rz(-theta)$)
- freeprecess(T,T1,T2,df) : calcul des opérateurs de relaxation en fonction du T1 et T2, de l'off-resonance de l'isochromat df durant un délai T

Validez que les fonctions fonctionnes bien (affichage d'une flèche avant et après rotation par exemple)
"

# ╔═╡ 1ba334ec-801e-4534-8856-81523a222b22
begin
	#--------------------------------------------------------
	#  By convention all rotations are left-handed
	#--------------------------------------------------------
	
	""" xrot(angle = 0., in_degs = False)
	Returns 4x4 matrix for left-handed rotation about x
	"""
	function xrot(angle = 0., in_degs = false)
	    if in_degs
	        angle = angle*π/180
		end
	    c = cos(angle)    
	    s = sin(angle)    
	    M = [[1,0,0,0] [0,c,-s,0] [0,s,c,0] [0,0.0,0.0,1]]'
	end
	
	""" yrot(angle = 0., in_degs = false)
	Returns 4x4 matrix for left-handed rotation about y
	"""
	function yrot(angle = 0., in_degs = false)
	    if in_degs
	        angle = angle*π/180
		end
	    c = cos(angle)    
	    s = sin(angle)    
	    M = [[c,0,s,0] [0,1,0,0] [-s,0,c,0] [0,0.0,0.0,1]]'
	end
	
	""" zrot(angle = 0., in_degs = false)
	Returns 4x4 matrix for left-handed rotation about y
	"""
	function zrot(angle = 0., in_degs = false)
	    if in_degs
	        angle = angle*π/180
		end
	    c = cos(angle)    
	    s = sin(angle)    
	    M = [[c,-s,0,0] [s,c,0,0] [0,0,1,0] [0,0.0,0.0,1]]'
	end
	
	""" freeprecess(t,T1 = 2000, T2 = 100, df = 0)
	
	Returns A function as a 4x4 matrix for freeprecession
	"""
	function freeprecess(t,T1 = 2000, T2 = 100, df = 0)
	    T1 = T1 * 1.
	    T2 = T2 * 1.
	    t = t * 1.
	    
	    phi = 2*pi*df*t/1000
	    E1 = exp(-t/T1)
	    E2 = exp(-t/T2)
	    A = diagm([E2, E2, E1,1])
		A[3,4] = 1-E1
		A = zrot(phi,false)*A
	    
		return A
	end
end;

# ╔═╡ 943d3397-b6fe-4ecf-ba21-a862d1a718b8
md"### Prepare some plot functions"

# ╔═╡ 863e9c02-a410-4135-8172-7c21605c58f0
""" plot_magnetization(M::Matrix{<:Real};title="")

M::Matrix{<:Real} : Magnetization vector of size 3(or 4)xN where N is the number of steps

Keywords :
- title of plots
"""
function plot_magnetization(M::Matrix{<:Real};title="",xlabel = "Time [ms]")
	f=Figure()
	ax = Axis(f[1,1],title=title)
	lines!(ax,M[1,:],label = "Mx")
	lines!(ax,M[2,:],label = "My")
	lines!(ax,M[3,:],label = "Mz")
	ax.xlabel= xlabel
	ax.ylabel="Mangnetization [ms]"
	axislegend()
	return f,ax
end

# ╔═╡ c080b7e5-4a24-46ea-bfa4-ab86a144d5f1
md"
## 2. Signal de précession libre

Simuler 2 signaux de precession libre (correspondant à l'évolution du signal après une impulsion radiofréquence)  avec une valeur d'off-resonance de 0 ou 10 Hz :

- angle = $$\frac{\pi}{2}$$ correspondant à une aimantation de départ selon l'axe 

$$\vec{M} = \begin{bmatrix}
  0\\ 
  1\\
  0 \\
  0
  \end{bmatrix}$$  

- Pas de simulation : 1 ms
- Durée de simulation : 1000 ms
- Off-resonance : 0 ou 10 Hz
- T1 = 600 ms
- T2 = 100 ms

**Expliquez la différence entre les 2 simulations ?**
"

# ╔═╡ a60b355b-d6f5-4e75-9b6b-99ef41b01759
function simuFID(df,dT,T,N::Int,T1,T2)
    ## Simulation with offresonance 
    A = freeprecess(dT,T1, T2, df)
    # Simulate the decay
	M = zeros(Float64,4,N)
    # initiliaze aimantation along x (like flip angle 90°)
    M[:,1] = [0,1,0,1]
    # propagate the relaxation along the time
    for k in range(1, N-1)
      M[:,k+1] = A * M[:,k]
	end
	return M
end

# ╔═╡ d0c30d9a-b6ea-43a8-99d1-f73752b9cc7c
begin
	# definition of parameters
	dT = 1
	T = 1000
	N = Int(ceil(T/dT))+1
	T1 = 600
	T2 = 100
end;

# ╔═╡ 79073dc4-0c13-48fc-94b7-544d4dd62059
@bind df_slider PlutoUI.Slider(0:0.5:50, default=0)

# ╔═╡ 5d15c345-ffef-4b4a-a628-b895096a59e4
 M = simuFID(df_slider,dT,T,N,T1,T2)

# ╔═╡ 3f9fcad3-f4fa-4d96-b41c-4437f2e2c65d
md"### Le slider controle la valeur d'off-resonance : $df_slider Hz"

# ╔═╡ e4bd5762-dbe7-47d3-83c0-afdf86efa41e
begin
	f=Figure()
	ax = Axis(f[1,1])
	lines!(ax,M[1,:],label = "Mx")
	lines!(ax,M[2,:],label = "My")
	lines!(ax,M[3,:],label = "Mz")
	ax.xlabel="Time [ms]"
	ax.ylabel="Mangnetization [ms]"
	ax.title="Evolution de l'aimantation (df = $df_slider)"
	
	axislegend()
	f
end

# ╔═╡ 3736e07f-8994-44e1-b27e-d0215a0b0e00
md"
## 3. Echo de gradient

![Gradient_Echo](img/Gradient_Echo.png)

Simuler le signal d'une séquence écho de gradient avec les paramètres ci-dessous :

```df = 0;	% Hz off-resonance.
df = 0;     % Hz
T1 = 600;	  % ms
T2 = 100;	  % ms
TE = 1;     % ms
TR = 500;   % ms
alpha = 60; % deg
```

* Quelle est la valeur de l'aimantation au temps TE après la 1ère impulsion de radiofréquence et après la seconde ? Expliquez la différence.
"

# ╔═╡ 867a2c0f-57f0-4121-bb04-647b15a894a7
begin
	df_2 = 0	# Hz off-resonance.
	T1_2 = 1400	# ms.
	T2_2 = 200	# ms.
	TE_2 = 1    #ms
	TR_2 = 100     #ms
	alpha_2 = pi/3      #deg
	Nex_2 = 50	# 20 excitations.
end;

# ╔═╡ 11a8d83e-f9ea-43d6-8c94-9037ee247510
begin
	# define usefull matrix
	Rflip = yrot(alpha_2)
	Atr = freeprecess(TR_2,T1_2, T2_2, df_2)
	Ate = freeprecess(TE_2,T1_2, T2_2, df_2)
	
	M0 = [0,0,1,1] # initialize along Z
	M_2 = Rflip * M0 # magnetization after flip
	M_2 = Ate * M_2 # relaxation to TE
	
	Mxy = sqrt(M_2[1]^2+M_2[2]^2)

	# after the second excitation (start back from M0)
	M_2 = Rflip * M0 # magnetization after flip
	M_2 = Atr * M_2 # relaxation for 1 TR
	M_2 = Rflip * M_2 # magnetization after flip
	M_2 = Ate * M_2 # relaxation to TE
	Mxy_2 = sqrt(M_2[1]^2+M_2[2]^2)
end

# ╔═╡ b2fdf1be-bf3f-4d20-aab3-b8066f3ebfa4
md"Signal for TR n°1 at TE  = $Mxy \
Signal for TR n°2 at TE  = $(Mxy_2)
"

# ╔═╡ 298dc120-a6f8-4a6b-969c-2e0d0c5ed784
md"* Simuler le signal pendant 10 TR et affichez l'évolution des aimantations selon x/y et z avec un pas dT = $dT ms. Faites de même avec seulement au temps TE."

# ╔═╡ d4fc5ffd-1676-4811-ab3d-f3534343a05f
function simuEchoDeGradient(df,alpha,TE,TR,dT,NEX,T1,T2)
    Ntr = round(TR/dT)
    M = zeros(Float32,4, floor(Int32,Ntr*NEX)+1) # store magnetization along the time
    Mte = zeros(Float32,4, NEX) # store magnetization for each TE

    At = freeprecess(dT,T1, T2, df)
    Ate = freeprecess(TE,T1, T2, df)

    Rflip = yrot(alpha)
    # initiliaze aimantation along z
    M[:,1] = [0,0,1,1]

    Mcount=1
    for n in range(1,NEX)
        M[:,Mcount] = Rflip * M[:,Mcount]
        Mte[:,n] = Ate * M[:,Mcount]

        for k in range(0,Ntr-1)
            Mcount = Mcount + 1
            M[:,Mcount] = At * M[:,Mcount-1]
		end
	end
    return M, Mte
end

# ╔═╡ 4f4b9886-cd8d-4686-8998-33c86a852b15
M2, Mte_2 = simuEchoDeGradient(df_2,alpha_2,TE_2,TR_2,dT,Nex_2,T1_2,T2_2)

# ╔═╡ d28f9c56-2365-4d82-9688-8837dea0cbd1
begin
	f2,_ = plot_magnetization(M2,title = "Evolution de l'aimantation au cours de 10 TR")
	f2
end

# ╔═╡ a137fb06-7333-4fd1-bf3e-05d3ce50e2f7
begin
	f3,_ = plot_magnetization(Mte_2,title = "Evolution de l'aimantation xy pour chaque TE",xlabel="Nombre d'excitation")
	f3
end

# ╔═╡ 00fe8d73-6db2-4a23-ac9a-c79858c68f13
md"* Combien de TR sont nécessaire pour atteindre un état d'équilibre ( **steady-state**) que l'on définira ici par une modification du signal inférieur à 1%."

# ╔═╡ 64e30d09-1fd5-496d-bc02-15939c859e63
begin
	for n = 2:size(Mte_2,2)
		if abs((Mte_2[1,n] - Mte_2[1,n-1])/Mte_2[1,n-1]) < 0.01
			display("The signal reach the steady states (< 1%) after $n excitations")
			break
		end
	end
end

# ╔═╡ 3fc59af9-9d32-49d7-9844-0c8bb4555e72
md"* Ecrire une fonction permettant de calculer le steady-state"

# ╔═╡ 03bc52e7-2df0-40c9-b2f2-5ec134c41927
begin
	function sssignal(flip,T1,T2,TE,TR,df)
	    # Steady-state directement au TE
	
	    # 	Starting at TE, M=M1
	    #	At TR, M=M2, and M2=Atetr*M1.
	    # 	At TE, M=M3, and M3=Ate*Rflip*M2.
	    #			M3=Ate*Rflip*Atetr*M1.
	    #	But M3=M1=Mss at steady state:
	    #   (I - Ate*Rflip*Atetr) * Mss = 0

	    Ate = freeprecess(TE,T1,T2,df)
	    Atetr = freeprecess(TR-TE,T1,T2,df)
	
		Aeq = Ate * Rflip * Atetr
	
	  	Mss = pinv((I-Aeq)[1:3,1:3])*(Aeq)[1:3,4]
	    Msig = complex(Mss[1],Mss[2])
	
	    return Msig, Mss
	end
	MsigTE,MssTE = sssignal(alpha_2,T1_2,T2_2,TE_2,TR_2,0)
	MsigTR,MssTR = sssignal(alpha_2,T1_2,T2_2,TR_2,TR_2,0)
end

# ╔═╡ 0013fdbb-37ea-4ea0-b8cb-1e143a764169
md"Le signal de steady state au temps TE au temps TE pour $(round(abs(MsigTE),digits=4))

Au *TR* : 
- Mx = $(round(MssTR[1],digits = 4))
- Mz = $(round(MssTR[3],digits = 4))
"

# ╔═╡ a4f4571d-8006-43e9-b6b8-8a4fe83aac7b
md"### Cas de la séquence FLASH

Dans le cas de la séquence FLASH les valeurs d'aimantations Mx et My sont égale à 0 à la fin du TR. L'équation analytique du signal est simplifiée et devient :

$$M_{ss} = \frac{M_0(1 - E_1)}{1 - E_1 \cos(\Theta)}$$
Qui donne  Mss = $((1-exp(-TR_2/T1_2))/(1-exp(-TR_2/T1_2)*cos(alpha_2)))

* Calculer la même chose pour une séquence spoilée à la fin du TR (FLASH). "

# ╔═╡ c53c055d-91de-41a3-8939-23a2e3c33ed0
begin
	function flashsignal(flip,T1,T2,TE,TR,df)
	    # Steady-state directement au TE
	
	    # 	Starting at TE, M=M1
	    #	At TR, M=M2, and M2=Atetr*M1.
	    # 	At TE, M=M3, and M3=Ate*Rflip*M2.
	    #			M3=Ate*Rflip*Atetr*M1.
	    #	But M3=M1=Mss at steady state:
	    #   (I - Ate*Rflip*Atetr) * Mss = 0

	    Ate = freeprecess(TE,T1,T2,df)
	    Atetr = freeprecess(TR-TE,T1,T2,df)

		spoiler = [[0 0 0 0];[0 0 0 0];[0 0 1 0];[0 0 0 1]]
		Aeq = Ate * Rflip * spoiler*Atetr
	
	  	Mss = pinv((I-Aeq)[1:3,1:3])*(Aeq)[1:3,4]
	    Msig = complex(Mss[1],Mss[2])
	
	    return Msig, Mss
	end
	MsigTR2,MssTR2 = flashsignal(alpha_2,T1_2,T2_2,TR_2,TR_2,0)
	MsigTE2,MssTE2 = flashsignal(alpha_2,T1_2,T2_2,TE_2,TR_2,0)
end

# ╔═╡ afd94a4c-d097-4a33-a651-2dd0dcc80e68
md"Le signal de steady state au temps TE au temps TE pour $(round(abs(MsigTE2),digits=4))

Au *TR* : 
- Mx = $(round(MssTR2[1],digits = 4))
- Mz = $(round(MssTR2[3],digits = 4))
"

# ╔═╡ f8547a0f-1b1b-4cd7-b0ee-51c58ff52603
md"# 4. Echo de spin

```julia
T1 = 600  #ms
T2 = 100  #ms
df = 10   #Hz
TR = 500  #ms
TE = 50   #ms
dT = 1    #ms 
```

Simuler le signal au cours d'un TR et affichez l'aimantation Mx, My et Mz. 
"

# ╔═╡ c4386a6f-b930-4751-ab59-026f1406d415
## Simulation Echo de spin sur un TR
function SpinEchoEvolution(;dT = 1,		# 1ms delta-time.
	T = 1000,	# total duration
	df = 10,	# Hz off-resonance.
	T1 = 600,	# ms.
	T2 = 100,	# ms.
	TE = 50,	# ms.
	TR = 500)	# ms.
	
	N1 = round(Int32,TE/2/dT)
	N2 = round(Int32,(TR-TE/2)/dT)
	
	A = freeprecess(dT,T1,T2,df)
	M = zeros(Float32,4, N1+N2+1)
	M[:,1]=[0,0,1,1]
	
	Rflip = yrot(pi/2)
	Rrefoc = xrot(pi)
	
	M[:,2]= A * Rflip * M[:,1]
	
	for k in range(3, N1+1)
	    M[:,k] = A * M[:,k-1]
	end
		
	M[:,N1+2]= A * Rrefoc * M[:,N1+1]   
	
	for k in range(2, N2)
	    M[:,k + N1 + 1] = A * M[:,k + N1]
	end
	return M
end

# ╔═╡ 196b0950-797c-44cf-8dcb-f72b8b62559c
begin
	f4,ax4 = plot_magnetization(SpinEchoEvolution();title="Aimantation of a spin echo")
	vlines!(ax4,50,color=:black)
	ax4.xticks=([0,50,collect(100:100:1000)...],["0","TE",string.(100:100:1000)...])
	f4
end

# ╔═╡ 9b9a4ba3-f0e3-4721-b5d5-fe8137dfaba5
md"Simulez 10 isochromates avec des valeurs différentes d'off-resonance (df) répartie de -50 à 50 Hz.

Affichez sur le même graphique: La valeur absolue d'un isochromate dans le plan transverse $M_{Sig}=\sqrt{x^2+y^2}$ et la valeur absolue de la somme des isochromates.

Et dans un second graphique, la phase de tout les isochromates au cours du temps
"


# ╔═╡ 24761ede-0303-424f-92ad-6df5add869bb
function simuEchoSpinMultiDf(;nDf::Int64 = 10)
	f=Figure()
	ax = Axis(f[1,1],title="evo")
	ax.xlabel="Time [ms]"
	ax.ylabel="Mangnetization [ms]"

	M = SpinEchoEvolution(df=0)
	M = zeros(eltype(M),size(M)...,nDf)

	randFreq = 50*rand(nDf).-50/2
	for (i,df) in enumerate(randFreq)
		M[:,:,i] = SpinEchoEvolution(df=df)
	end
	lines!(ax,sqrt.(M[1,:,1].^2+M[2,:,1].^2),label = "Mxy")
	
	Mave = mean(M,dims=3)[:,:,1]
	lines!(ax,sqrt.(Mave[1,:].^2+Mave[2,:].^2),label = "Mxy of sum")
	axislegend()

	vlines!(ax,50,color=:black)
	ax.xticks=([0,50,collect(100:100:1000)...],["0","TE",string.(100:100:1000)...])
	
	# plot phase
	ax2 = Axis(f[2,1],title="Phase evolution of all the isochromate")
	ax2.xlabel="Time [ms]"
	ax2.ylabel="Phase [rad]"
	
	for (i,df) in enumerate(randFreq)
		lines!(ax2,angle.(M[1,:,i]+im*M[2,:,i]))
	end
	vlines!(ax2,50,color=:black)
	ax2.xticks=([0,50,collect(100:100:1000)...],["0","TE",string.(100:100:1000)...])
	
	f
end


# ╔═╡ 9b8b5da9-388a-4cf2-b208-712bde4f650b
simuEchoSpinMultiDf(nDf = 10)

# ╔═╡ 00fd968e-b63e-4780-b923-91a94e1258d7
md"Ecrivez une fonction pour calculer l'état d'équilibre de la séquence écho de spin"

# ╔═╡ f207da71-5f44-4691-9671-d9673133e267
"""sesignal(;df = 0,	# Hz off-resonance.
T1 = 600,	# ms.
T2 = 100,	# ms.
TE = 50,	# ms.
TR = 1000)	# ms.
"""
function sesignal(;df = 0,	# Hz off-resonance.
T1 = 600,	# ms.
T2 = 100,	# ms.
TE = 50,	# ms.
TR = 1000)	# ms.
	Rflip = yrot(pi/2) # Rotation from excitation pulse (90)
	Rrefoc = yrot(pi) # Rotation from refocusing pulse (usually 180)
	Atr = freeprecess(TR-TE,T1,T2,df) # Propagation TE to TR
	Ate2,Bte2 = freeprecess(TE/2,T1,T2,df) # Propagation 0 to TE/2 (same as TE/2 to TE)


	spoiler = [[0 0 0 0];[0 0 0 0];[0 0 1 0];[0 0 0 1]]
	Aeq = Ate2 * Rrefoc * Ate2 * Rflip * spoiler * Atr

	Mss = pinv((I-Aeq)[1:3,1:3])*(Aeq)[1:3,4]
	Msig = complex(Mss[1],Mss[2])

	return Msig, Mss
end

# ╔═╡ 010f2f77-5888-424d-8c3e-bf6e9f1b4263
begin
	Msig,Mss = sesignal(df = 0,T1 = 600,T2 = 100,TE = 50,TR = 1000)
	abs(Msig)
end

# ╔═╡ 732ca04e-c7f5-423f-a1da-3ce576bcd514
md"# Application : Optimisation du contrast sur un cerveau"

# ╔═╡ 3f93ed2a-be66-402c-8d98-7ff14c3b67e7
md"Nous allons visualiser l'effet des choix de paramètre sur une séquence. Le jeu de donnée que nous allons charger est un cerveau dont les valeurs sont comprises entre 0 et 4:
-0 : vide
- 1 : WM
- 2 : GM
- 3 : LCR"

# ╔═╡ 2f44e437-6247-4eb7-b2a6-b621256ad11c
begin
	mask=niread("mask.nii").raw
	heatmap(mask,colormap=:greys)
end

# ╔═╡ 2f426a74-52bf-45a9-aa5b-96268dc495f9
begin
	mask2 = zeros(Int32,size(mask))
	mask2[mask .== 1] .=1
	heatmap(mask2,colormap=:greys)
end

# ╔═╡ e10aa737-92e8-4db5-b426-dbba71e36ebc
md"Paramètres de simulation"

# ╔═╡ dff344d0-b0a9-41e7-8d8a-b684afe895c7
begin
	# tissues params
	T1_mask=[1400,800,3000]
	T2_mask=[40,70,2000]
end

# ╔═╡ cc607d22-408e-4c4d-ac93-e95ee52bce85
@bind TR_slider PlutoUI.Slider(0:10:3000, default=1000)

# ╔═╡ 7631ad90-6055-44df-bb68-9154f49bdeec
md"TR = $TR_slider ms"

# ╔═╡ 577e26b1-bc7b-40e0-b04b-85f46071c421
@bind TE_slider PlutoUI.Slider(0:1:200, default=80)

# ╔═╡ fc1f3f96-96ac-45d3-8e60-46e9dded075f
begin
	brain = zeros(Float32,size(mask))
	Mss_se = Float32[]
	for i =1:3
		Msig,Mss = sesignal(df = 0,T1 = T1_mask[i],T2=T2_mask[i],TE = TE_slider,TR = TR_slider)
		brain[mask .== i] .= abs(Msig)

		push!(Mss_se,abs(Msig))
	end
	CNR = abs(Mss_se[1]-Mss_se[2])/sqrt(TR_slider) # contrast normalized by time
end

# ╔═╡ 2c02540c-4af9-44ae-b4d6-2e68a09cfd6d
md"TE = $TE_slider ms"

# ╔═╡ deb9e407-b0fd-4a66-a869-a4de97471f4a
begin
	f6,ax6,h6 = heatmap(brain,colormap=:greys)
	ax6.title="CNR = $CNR"
	hidedecorations!(ax6)
	f6
end

# ╔═╡ 5848223b-262d-4a23-b123-84616e074fb1
md"## CNR optimum"

# ╔═╡ 74b1c8f4-6ce5-46c6-80e4-c7765088f992
begin
	CNR_mat = zeros(Float32,200,500)
	range_TR = LinRange(10,5000,500)
	range_TE = LinRange(5,200,200)
	
	for (i,TR) in enumerate(range_TR)
		for (j,TE) in enumerate(range_TE)

			if TE <= TR
			Msig,_ = sesignal(df = 0,T1 = T1_mask[1],T2=T2_mask[1],TE = TE,TR = TR)
			Msig2,_ = sesignal(df = 0,T1 = T1_mask[2],T2=T2_mask[2],TE = TE,TR = TR)
	
			CNR_mat[j,i] = abs(Msig-Msig2)/sqrt(TR)
			else
				CNR_mat[j,i] = NaN
			end
		end
	end
end

# ╔═╡ ea594a6f-2092-45d2-be86-a8e85167f376
begin
	CNR_mat2 = copy(CNR_mat)
	CNR_mat2[isnan.(CNR_mat)] .= 0
	val,index = findmax(CNR_mat2)
	
	f7,ax7,h7 = heatmap(CNR_mat)
	ax7.title="Maximum CNR = $val with TE = $(range_TE[index[1]]) and TR = $(range_TE[index[2]])"
	ax7.xlabel= "TE [ms]"
	ax7.ylabel="TR [ms]"

	scatter!(ax7,index[1],index[2],color=:red,marker = :cross,markersize = 30)
	f7
	
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
NIfTI = "a3a9e032-41b5-5fc4-967a-a6b7a19844d3"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CairoMakie = "~0.10.2"
NIfTI = "~0.4.0"
PlutoUI = "~0.7.49"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "bf63cda348d3b58e21474de664483d09dd20adb0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0310e08cb19f5da31d08341c6120c047598f5b9c"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.5.0"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "SnoopPrecompile"]
git-tree-sha1 = "abb7df708fe1335367518659989627100a61f3f0"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.10.2"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "844b061c104c408b24537482469400af6075aae4"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.5"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "74911ad88921455c6afcad1eefa12bd7b1724631"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.80"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "d3ba08ab64bdfd27234d3f61956c966266757fe6"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.7"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "57f7cde02d7a53c9d1d28443b9f11ac5fbe7ebc9"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.3"

[[deps.GZip]]
deps = ["Libdl"]
git-tree-sha1 = "039be665faf0b8ae36e089cd694233f5dee3f7d6"
uuid = "92fee26a-97fe-5a0c-ad85-20a5f3185b63"
version = "0.5.1"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "e07a1b98ed72e3cdd02c6ceaab94b8a606faca40"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.2.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "fe9aea4ed3ec6afdfbeb5a4f39a2208909b162a6"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.5"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "678d136003ed5bceaab05cf64519e3f956ffa4ba"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "106b6aa272f294ba47e96bd3acbabdc0407b5c60"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.2"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "680e733c3a0a9cea9e935c8c2184aea6a63fa0b5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.21"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "MiniQhull", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Setfield", "Showoff", "SignedDistanceFields", "SnoopPrecompile", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "274fa9c60a10b98ab8521886eb4fe22d257dca65"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.2"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "2c3fc86d52dfbada1a2e5e150e50f06c30ef149c"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.2"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test", "UnicodeFun"]
git-tree-sha1 = "f04120d9adf4f49be242db0b905bea0be32198d1"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.4"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.MiniQhull]]
deps = ["QhullMiniWrapper_jll"]
git-tree-sha1 = "9dc837d180ee49eeb7c8b77bb1c860452634b0d1"
uuid = "978d7f02-9e05-4691-894f-ae31a51d76ca"
version = "0.4.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NIfTI]]
deps = ["Base64", "GZip", "Mmap", "Test"]
git-tree-sha1 = "2151282446bb161941bff1b330a3212fb132e6b1"
uuid = "a3a9e032-41b5-5fc4-967a-a6b7a19844d3"
version = "0.4.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "5ae7ca23e13855b3aba94550f26146c01d259267"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "f809158b27eba0c18c269cf2a2be6ed751d3e81d"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.17"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "18f84637e00b72ba6769034a4b50d79ee40c84a9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.5"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "c95373e73290cf50a8a22c3375e4625ded5c5280"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QhullMiniWrapper_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Qhull_jll"]
git-tree-sha1 = "607cf73c03f8a9f83b36db0b86a3a9c14179621f"
uuid = "460c41e3-6112-5d7f-b78c-b6823adb3f2d"
version = "1.0.0+1"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "238dd7e2cc577281976b9681702174850f8d4cbc"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1001+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "786efa36b7eff813723c4849c90456609cf06661"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "8b20084a97b004588125caebf418d8cab9e393d1"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.4"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StableHashTraits]]
deps = ["CRC32c", "Compat", "Dates", "SHA", "Tables", "TupleTools", "UUIDs"]
git-tree-sha1 = "0b8b801b8f03a329a4e86b44c5e8a7d7f4fe10a3"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "0.3.1"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "129703d62117c374c4f2db6d13a027741c46eafd"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.13"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "b03a3b745aa49b566f128977a7dd1be8711c5e71"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.14"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "7e6b0e3e571be0b4dd4d2a9a3a83b65c04351ccc"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.3"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

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

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

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
"""

# ╔═╡ Cell order:
# ╟─11bf61dc-aa97-48b6-ab67-9bc0861120bc
# ╠═2e4cbd0e-a92b-11ed-022c-eda67b72483f
# ╠═1ba334ec-801e-4534-8856-81523a222b22
# ╟─943d3397-b6fe-4ecf-ba21-a862d1a718b8
# ╠═863e9c02-a410-4135-8172-7c21605c58f0
# ╠═c080b7e5-4a24-46ea-bfa4-ab86a144d5f1
# ╠═a60b355b-d6f5-4e75-9b6b-99ef41b01759
# ╠═d0c30d9a-b6ea-43a8-99d1-f73752b9cc7c
# ╠═5d15c345-ffef-4b4a-a628-b895096a59e4
# ╟─3f9fcad3-f4fa-4d96-b41c-4437f2e2c65d
# ╟─79073dc4-0c13-48fc-94b7-544d4dd62059
# ╠═e4bd5762-dbe7-47d3-83c0-afdf86efa41e
# ╟─3736e07f-8994-44e1-b27e-d0215a0b0e00
# ╠═867a2c0f-57f0-4121-bb04-647b15a894a7
# ╠═11a8d83e-f9ea-43d6-8c94-9037ee247510
# ╟─b2fdf1be-bf3f-4d20-aab3-b8066f3ebfa4
# ╟─298dc120-a6f8-4a6b-969c-2e0d0c5ed784
# ╠═d4fc5ffd-1676-4811-ab3d-f3534343a05f
# ╠═4f4b9886-cd8d-4686-8998-33c86a852b15
# ╠═d28f9c56-2365-4d82-9688-8837dea0cbd1
# ╠═a137fb06-7333-4fd1-bf3e-05d3ce50e2f7
# ╟─00fe8d73-6db2-4a23-ac9a-c79858c68f13
# ╠═64e30d09-1fd5-496d-bc02-15939c859e63
# ╟─3fc59af9-9d32-49d7-9844-0c8bb4555e72
# ╠═03bc52e7-2df0-40c9-b2f2-5ec134c41927
# ╟─0013fdbb-37ea-4ea0-b8cb-1e143a764169
# ╠═a4f4571d-8006-43e9-b6b8-8a4fe83aac7b
# ╠═c53c055d-91de-41a3-8939-23a2e3c33ed0
# ╟─afd94a4c-d097-4a33-a651-2dd0dcc80e68
# ╟─f8547a0f-1b1b-4cd7-b0ee-51c58ff52603
# ╠═c4386a6f-b930-4751-ab59-026f1406d415
# ╠═196b0950-797c-44cf-8dcb-f72b8b62559c
# ╟─9b9a4ba3-f0e3-4721-b5d5-fe8137dfaba5
# ╠═24761ede-0303-424f-92ad-6df5add869bb
# ╠═9b8b5da9-388a-4cf2-b208-712bde4f650b
# ╟─00fd968e-b63e-4780-b923-91a94e1258d7
# ╠═f207da71-5f44-4691-9671-d9673133e267
# ╠═010f2f77-5888-424d-8c3e-bf6e9f1b4263
# ╟─732ca04e-c7f5-423f-a1da-3ce576bcd514
# ╟─3f93ed2a-be66-402c-8d98-7ff14c3b67e7
# ╠═2f44e437-6247-4eb7-b2a6-b621256ad11c
# ╠═2f426a74-52bf-45a9-aa5b-96268dc495f9
# ╠═e10aa737-92e8-4db5-b426-dbba71e36ebc
# ╠═dff344d0-b0a9-41e7-8d8a-b684afe895c7
# ╠═fc1f3f96-96ac-45d3-8e60-46e9dded075f
# ╟─7631ad90-6055-44df-bb68-9154f49bdeec
# ╠═cc607d22-408e-4c4d-ac93-e95ee52bce85
# ╟─2c02540c-4af9-44ae-b4d6-2e68a09cfd6d
# ╠═577e26b1-bc7b-40e0-b04b-85f46071c421
# ╟─deb9e407-b0fd-4a66-a869-a4de97471f4a
# ╟─5848223b-262d-4a23-b123-84616e074fb1
# ╠═74b1c8f4-6ce5-46c6-80e4-c7765088f992
# ╠═ea594a6f-2092-45d2-be86-a8e85167f376
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
