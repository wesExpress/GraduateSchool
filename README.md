# GraduateSchool
Various code that I have written over my graduate career to analyze a wide range of galaxy data.


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
This python routine reads in two on-sky galaxy image as text files and computes three different colors. In astronomy, we observe objects in different wavelength bands, as certain stars, galaxies etc. give off light differently depending on the wavelength. We simply define a color as the difference between two photometric bands. Here, we are using two images, one through a blue filter (Johnson B) and one through a red filter (Johnson I). 

The three different colors we measure are all *B-I*. However, each take different components of the galaxy into account. 
-**bar color**: defined within the bar region
-**disk color**: defined outside the bar region
-**area color**: defined as the entire galaxy out to a certain radius
As bars and disks have different stellar populations, we can learn a lot about a galaxy by looking at the colors in each of these regions.

## deproject.py
Similar to bar_len.py above, this python routine analyzes elliptical isophotes. Here, we simply obtain the position angle and inclination of the galaxy.

## imread.py
A simple python routine that reads in a \*.fits file and converts it to a text file.

## phase.f90
Another backbone of my thesis project.

This fortran routine requires azimuthal light profiles output by **az_prof.f90**. Here, we decompose the azimuthal light profiles via a Fourier transform in order to learn more about the bar. As a bar is at a constant azimuthal angle, it shows up twice in 0 to 2pi azimuthal decomposition. Therefore, a Fourier analysis will show that a bar has very strong even modes of a Fourier series. This is an additional means of proving a bar is in a galaxy. In addition, there is also a means of measuring a bar length via the Fourier Amplitudes. Finally, we use this routine to measure an additional feature of a galaxy bar: the corotation radius. This is simply the radius where orbits in a galaxy are equal to the pattern speed of the bar. This is important for determining things like the relative bar pattern speed, and is one of the three main bar properties measured in astronomy.
