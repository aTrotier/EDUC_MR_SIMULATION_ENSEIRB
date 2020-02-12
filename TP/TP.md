
<center><h1>TP IRM SIMULATION</h1>
<b>Aurélien Trotier</b> : Février 2020
</center>


Le signal IRM peut-être représenté grâce à l'équation de Bloch :

$$\frac{d\vec{M}}{dt} = \gamma \vec{M} \times \vec{B}_{ext} + \frac{1}{T_1} (M_0 - M_z)\vec{z} - \frac{1}{T_2} \vec{M}_\perp$$

<h3>1. Préparation des sous-fonctions</h3>
Préparez les fonctions nécessaires à la simulation du signal IRM :
- xrot(angle), yrot(angle) et zrot(angle) : qui retourne une matrice de rotation à partir d'un angle (en radian)
- throt(phi,theta) : si l'angle de rotation est transverse (différent de x/y) (utilisez les fonctions précédentes : $Rth(phi,theta)=Rz(theta)*Rx(phi)*Rz(-theta)$)
- freeprecess(T,T1,T2,df) : calcul des opérateurs de relaxation en fonction du T1 et T2, de l'off-resonance de l'isochromat df durant un délai T

Validez que les fonctions fonctionnes bien (affichage d'une flèche avant et après rotation par exemple)

<h3>2. Signal de précession libre</h3>
Simuler 2 signaux de precession libre (correspondant à l'évolution du signal après une impulsion radiofréquence)  avec une valeur d'off-resonance de 0 ou 10 Hz :

- angle = $$\frac{\pi}{2}$$ correspondant à une aimantation de départ selon l'axe $$\vec{M} = \begin{bmatrix}
   0\\ 
   1\\
  0
  \end{bmatrix}$$  
- Pas de simulation : 1 ms
- Durée de simulation : 1000 ms
- Off-resonance : 0 ou 10 Hz
- T1 = 600 ms
- T2 = 100 ms

<u>Expliquez la différence entre les 2 simulations ?</u>



<h3>3. Echo de gradient</h3>
![Gradient_Echo](/Users/aurelien/Documents/WORK/Cours/ENSEIRB/2020_ENSEIRB_MR_SIMULATION/TP/img/Gradient_Echo.png)

Simuler le signal d'une séquence écho de gradient avec les paramètres ci-dessous :

```matlab
df = 0;     % Hz
T1 = 600;	  % ms
T2 = 100;	  % ms
TE = 1;     % ms
TR = 500;   % ms
alpha = 60; % deg
```



* Quelle est la valeur de l'aimantation au temps TE après la 1ère impulsion de radiofréquence et après la seconde ? Expliquez la différence.

* Combien de TR sont nécessaire pour atteindre un état d'équilibre ( **steady-state**) que l'on définira ici par une modification du signal inférieur à 1%.
* Calculez la valeur exact de l'état d'équilibre au TR puis au TE. Ecrivez une fonction **[Msig,Mss]=sssignal(flip,T1,T2,TE,TR,dfreq) ** permettant de calculer le signal de steady-state au TE (Mss) et juste la valeur complexe dans le plan : $$Msig = Mss(1)+i*Mss(2)$$

Dans le cas d'une séquence **FLASH** (**F**ast **L**ow **A**ngle **SH**ot) l'aimantation steady-state longitudinal (selon l'axe z) est égale à $$M_{ss}= \frac{M_0(1-E_1)}{1-E_1 \cos \theta} = 0.7224$$ . Cette différence s'explique par le fait que l'aimantation est spoilée à la fin de chaque TR (RF et phase spoiling) donc que l'aimantation transversale (selon x/y) est égale à 0.

* Implémentez le *spoiling parfait*  et vérifiez que la valeur correspond à l'équation analytique ci-dessus. Ecrire la fonction **[Msig,Mss]=srsignal(flip,T1,T2,TE,TR,dfreq) ** correspondante.



<h3>4. Spin echo</h3>
<h5><u>Simulation</u></h5>


![Spin_Echo](/Users/aurelien/Documents/WORK/Cours/ENSEIRB/2020_ENSEIRB_MR_SIMULATION/TP/img/Spin_Echo.png)

Simuler le signal au cours d'un TR et affichez l'aimantation Mx, My et Mz. 

```matlab
T1 = 600  ms
T2 = 100  ms
df = 10   Hz
TR = 500  ms
TE = 50   ms
dT = 1    ms 
```



* Qu'observez vous au temps TE ? Comment être sur que c'est bien un écho de spin ?

Simulez 5 isochromates avec des valeurs différentes d'off-resonance (df) puis affichez  sur le même graphique: La valeur absolue d'un isochromate dans le plan transverse ($M_{Sig}=\sqrt{x^2+y^2}$) la valeur absolue de la somme des isochromates. 

* Expliquez pourquoi la courbe de la somme (divisée par le nombre d'isochromates) croise la courbe de la valeur absolue.

Refaites la même simulation plusieurs fois pour voir les différences et augmenter le nombre d'isochromates dans la simulation



<h5><u>Steady-state</u></h5>
* Calculez le steady state au TE de la séquence Spin-Echo

* Etendre ce calcul pour une séquence Turbo-spin echo (Echo Train Lenght = 8) [voir présentation pour propagation des équations].

  ![Turbo_Spin_Echo](/Users/aurelien/Documents/WORK/Cours/ENSEIRB/2020_ENSEIRB_MR_SIMULATION/TP/img/Turbo_Spin_Echo.png)

* Calculez le steady-state au 1er TE et comparez avec le steady-state de la séquence Spin-Echo. Vérifiez la valeur en mettant un Echo Train Length = 1. Expliquez ce résultat.




<h3>OPTIONAL : Gradient Spoiling et RF-Spoiling</h3>

Le spoiling consiste à utiliser des méthodes pour que l'aimantation transverse soit égale à 0 avant la prochaine impulsions de radio-fréquence. Il existe 3 méthodes de spoiling :

1. Spoiling naturel : TR >> T2 du tissus
2. Gradient spoiling : Dephasage des spins dans un voxel à l'aide d'un gradient
3. RF-Spoiling + Gradient spoiling: L'impulsion de radiofréquence utilisé pour chaque répétition est différent 



<h5><u>Question 1</u></h5>

Comment vérifier que le spoiling naturel donne le même résultat que le spoiling parfait que nous avons utilisé précédement ?

<h5><u>Question 2</u></h5>

Un spoiler gradient consiste à déphaser les spins dans un voxel avant la prochaine répétition. Généralement grâce à un fort gradient selon un axe. Pour le simuler il va falloir simuler plusieurs isochromats dans un voxel qui vont être déphasé différement par le gradient puis en faire la somme.

* Ecrire une fonction **Mss=gssignal(flip,T1,T2,TE,TR,dfreq,phi)** ou phi correspond à l'angle de dephasage à la fin du TR

* Vérifiez que la valeur est égale à Mss = [0.1248, 0.1129, 0.1965]' avec les paramètres

  ```matlab
  flip = pi/3;
  T1 = 600;
  T2 = 100;
  TE = 2;
  TR = 10;
  dfreq = 0;
  phi= pi/2;
  ```

* Ecrire une fonction **Mss=gresignal(flip,T1,T2,TE,TR,dfreq)** permettant de calculer l'aimantation moyenne dans un voxel pour 100 isochromates qui subissent un spoil de -2pi à 2pi. Vérifiez la valeur du signal avec les même paramètres : Mss=[0.1157, 0.0000, 0.1801]'
* Comparez le signal



<h5><u>Question 3</u></h5>

Le RF spoiling consiste exciter l'aimantation à chaque TR avec une phase différente (en utilisant **throt(flip,Rfph)**). Généralement on utilise : 

```matlab
inc = 117;
Rfph(1)=0;
for n = 2:5
    Rfph(n)=(n-1)*inc+Rfph(n-1)
end
```

 Donc les 5 premières excitations auront une phase de : 0/117/351/792/1170 °. 

Attention la phase de la dernière excitation utilisé est soustraite à la phase du signal recueilli.

* Simulez et affichez le signal en magnitude et en phase en fonction du nombre d'excitation pour les 100 premières excitations. Qu'observez-vous ?

```matlab
flip = pi/6;
T1 = 600;
T2 = 100;
TE = 2;
TR = 10;
dfreq = 0;
```

* Comparer les valeurs à l'état d'équilibre de la séquence RF-spoiling avec la séquence ayant un spoiling parfait.
* Vérifiez que les valeurs sont proches pour toutes les valeurs de flip angle entre la séquence rf-spoiling (spgrsignal), spoiling parfait (srsignal) et gradient-spoiling (gresignal). Pour cela écrire une fonction : ***Msig=spgrsignal(flip,T1,T2,TE,TR,dfreq,Nex,inc)*** avec le nombre d'exitation Nex = 100 et inc = 117. Et faites varier l'angle de 0 à 90 degrès. Puis affichez les 3 courbes.

<h3>OPTIONAL : BSSFP et banding artefacts</h3>
BSSFP est l'acronyme pour "Balanced Steady-State Free Precession". C'est une dans laquelle les gradients sont équilibré (dont la somme au cours du TR = 0)

<img src="/Users/aurelien/Documents/WORK/Cours/ENSEIRB/2020_ENSEIRB_MR_SIMULATION/TP/img/bssfp_sequence.png" alt="bssfp_sequence" style="zoom: 50%;" /><img src="/Users/aurelien/Documents/WORK/Cours/ENSEIRB/2020_ENSEIRB_MR_SIMULATION/TP/img/bssfp.png" alt="bssfp" style="zoom:50%;" />

Cette séquence est donc équivalente à la séquence d'écho de gradient simulé dans la partie 3. Elle a pour particularité d'avoir un fort signal qui est un compromis T1/T2 avec des TRs court ainsi que de généré des **artéfacts de bande noire**.

Pour le simuler il faut prendre un compte l'off-resonance qui peut exister dans un voxel. 



