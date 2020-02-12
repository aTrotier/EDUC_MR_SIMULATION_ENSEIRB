<center><h1>PROJET IRM SIMULATION</h1>
<b>Aurélien Trotier</b> : Février 2020
</center>



<h3>Objectif</h3>

Générer un jupyter-notebook permettant de générer une image de cerveau acquise avec des paramètres spécifique pour une séquence. 

Le langage utilisé sera **python** et les graphiques devront être interactif (voir cet [exemple](https://qmrlab.org/jekyll/2018/10/23/T1-mapping-inversion-recovery.html))

<h3>Ressource disponibles</h3>

- Carte T1 du cerveau
- Carte T2 du cerveau
- Pas de carte M0 (donc on considère que tous les voxels ont la même quantité d'aimantation initiale)
- Equation du signal

<h3>Séquence à simuler</h3>

<h5><u>Séquence Flash</u></h5>

La séquence Flash est un écho de gradient, pour pouvoir visuliser les effets des spoilers il faut simuler différents isochromat dans un voxel qui font se déphaser de manière différentes dans celui-ci car ils ne sont pas à la même position physique.

Echo de gradient avec différent spoiler :

* Choix : spoiler RF

* Choix spoiler parfait $M_{xy}=0$

* Intensité spoiler gradient (de 0 à 10 $\pi$  avec un pas de $\pi/100$; valeur initiale : $2 \pi$)

  

Paramètres à modifier : 

* TE
* TR
* Angle
* nombre d'isochromat à simuler

<span style='color:red'>ATTENTION : Les limites devront être prise en compte (exemple : TR minimum limité par le TE) </span>



<h5><u>Spin Echo</u></h5>

L'équation analytique de la séquence peut-être utilisé : 

![Eq_Spin_echo](/Users/aurelien/Documents/WORK/Cours/ENSEIRB/2020_ENSEIRB_MR_SIMULATION/TP/img/Eq_Spin_echo.png) 

Cette équation néglige le signal selon $M_{xy}$ soit parce que le TR >> T2 des tissues ou grâce à un gradient



Paramètres à modifier : 

* TE (de $\approx 0$ à 10 $\pi$  à 1000 ms avec un pas de 5 ms)
* TR (de $\approx 0$ à 10 $\pi$  à 10000 ms avec un pas de 10 ms)
* Angle excitation $\theta_1$ et refocalisation $\theta_2$ (de 0 à 270 avec un pas de 1 degrès)

<span style='color:red'>ATTENTION : Les limites devront être prise en compte (exemple : TR minimum limité par le TE </span>



<h3>Résultats</h3>

Le jupyter notebook devra m'être envoyé avant **le 2 avril ** à l'adresse : <span style='color:blue'><u>trotier@rmsb.u-bordeaux.fr</u></span>


