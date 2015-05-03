**NOT YET TESTED**. Use at your own risk.

# Alignment #

The alignment scheme is similar to the one used in EQMOD: the driver maintains a set of alignment points and builds a model using equations (5.4.1) to (5.4.4) of the [Taki's paper](http://www.geocities.jp/toshimi_taki/index.htm). Alignment points are by default added to the model using the sync mechanism. The driver uses constantly the model to transform true mount coordinates given by the encoders into aligned celestial coordinates, which are then displayed in a planetarium as kstars. Conversely celestial coordinates as given by the planetarium when performing a goto are also transformed into mount coordinates according to the model. These two transformations use three points among the set of alignment points to build a 3x3 matrix and its inverse, and just multiply the suitable matrix with input coordinates to produce their results. Given an input coordinate, these three points are chosen according to the Delaunay triangulation of the set of alignment points: the model determines which triangle the input point belongs to, and uses its three vertices to build the matrices. If there is no such triangle the model uses a nearest point algorithm. The Delaunay triangulation is performed using a 3D convex hull algorithm as alignment points are on the unit sphere (the origin is first added to the set of points, then all faces containing it are removed from the convex hull). The source code for the Delaunay triangulation comes from the book of Joseph O'Rourke [Computational Geometry in C](http://cs.smith.edu/~orourke/books/ftp.html). The use of a Delaunay triangulation is a difference with the alignment model of the EQMOD project, where a triangle of alignment points is selected  so that the distance of its center to the target is minimal.

Besides this alignment model, the standard syncs mechanism is still active: a standard synchronisation simply shifts the alignment model. They may be useful for external tools which rely on such mechanism. In this case it is advisable to deactivate the alignment model ('No Alignment' or 'Clear List').

Alignment points may be saved to and restored from a 'Data file': this is useful for permanent observatories. The data file is automatically loaded when you connect the mount, which means that if the file exists, you use a permanent observatory. This file has a simple and self-explanatory XML content, so you could edit it by hand.

## Steps for using N-Star alignment in Kstars ##
  1. The scope location and the UTC time should be **correctly** set in the driver: you could tell kstars to set these two values when connecting devices in its option tab.
  1. Connect the mount.
  1. Add alignement points:
    * Track to a star of your choice
    * Align it in a reticle using 'Move N/S' and 'Move 'E/W' of the motion tab
    * Use 'Sync' in kstars on the star you choosed to add an alignment point
    * At least three alignment points should be defined (otherwise nearest point is choosen)
    * You could empty the set at any time if you wish to restart

Subsequent Track/Slew commands will now be aligned, and the correction displayed in the Message window: you could check if it is safe and abort motion if not. Note that RA/DEC coords will now be slightly changing even when RA tracking, as your mount is not perfectly polar-aligned and the alignement model reflects this misalignment.

## Polar align helper ##

The driver can compute an estimation of the polar misalignment using the two star alignement method described in Rob Johnson [paper](http://www.whim.org/nebula/math/pdf/twostar.pdf). The driver also suggests in the message window a "Goto" target
which you could use to realign the second star using only the alt/az adjustments of the mount.
  1. Set "No Alignment" in "Align" Tab, and "Standard Syncs" in "Sync" Tab
  1. "Clear Delta" in "Sync" Tab
  1. Track, align and sync star1
  1. Track, align and sync star2
  1. An estimation of the alt/az misalignment is displayed in the "Sync" Tab
  1. In the message window a "Goto" target is also suggested:
    * track to that target :
      1. switch "On set" to "Track"
      1. copy/paste suggested RA/DEC coords to the EOD coords
      1. press the corresponding "Set" button: the mount should slew to that target
    * realign the second star using the **altitude and azimuth knobs** of the mount: do make any other moves to the mount (apart letting the RA tracking active)
    * when the second star is realigned, "Sync" (still standard) the mount to the second star: sync delta should be very low now and your mount reticle should be back near the second star in the planetarium
    * eventually clear the sync delta