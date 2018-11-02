# GraduateSchool
Various code that I have written over my graduate career.


## az_prof.f90
This fortran program is one the backbones of my thesis project. 

In order to study bars in galaxies, we decompose galaxy images into *azimuthal light profiles*. In simple terms, these are light profiles at a constant radius in a galaxy. Instead of starting in the center of a galaxy and examining the light profile moving outwards at say a single azimuthal angle, we obtain the light profile at all azimuthal angles for each radius. As galaxies are randomly oriented on the sky to our line of sight, our profiles will trace ellipses on the galaxy. In order to simply this, we must *de-project* the images by obtaing the position angle and inclination of the galaxy disk so that our azimuthal light profiles simply trace concentric circles about the galaxy center.

az_prof takes in the deprojected galaxy image as a text file (converted from a .fits file by a Python script) and all the information needed to contrsuct the azimuthal light profiles:
-angdiv: the number of azimuthal divisions (the smaller the better)
-raddiv: the number of radial divisions
-rstart: the radius to start constructing the profiles
-xcen,ycen: the pixel centers
-nx, ny: the pixel size of the image

As images are limited by their angular resolution, we cannot start contructing profiles at r = 0. Therefore, there is a balancing act of limiting the number of azimuthal and radial bins and starting too far out to lose information close in the galaxy.

The azimuthal light profiles are finally output at the end to be used in further analysis routines.

## az_prof_fit.py
This python routine reads in the output of az_prof.f90 to measure the bar length in galaxies, a very important feature in galaxy structure.

The general idea behind this routine is to think about a stellar bar in a galaxy. These are features in the centers of galaxies that remain at constant azimuthal angle for its entire length. Therefore, if we construct azimuthal light profiles of a barred galaxy, we can determine the bar length by locating the azimuthal centroid of a bar. If we construct the light profiles from 0 to 2pi, the bar will show up *twice*, as two 'humps' in the profile. Thus, all we must do is measure the azimuthal centroid of these humps for a number of radii and find where either **a)** the humps vanish and the light profile resembles that of a pure disk, or **b)** the humps move to a new azimuthal centroid for various reasons. The radius at which either of these occur is therefore the length of the bar.

## bar_len.py
This python routine measures the length of a galaxy bar by means of elliptical isophotes. These are simply ellipses fit to galaxy images that contain equal light throughout. By varrying the position angle and ellipticity of each ellipse, we can move out radially and obtain ellipticity and position angle profiles. With these, we can obtain the position angle and inclination of the galaxy disk, which allows us to deproject the image. 

To measure the bar length with these elliptical isophotes, this code analyzes the radial behaviour of both the ellipticity and position angle and uses two different methods of obtaining a bar length.

## colors.py
