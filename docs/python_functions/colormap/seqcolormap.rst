.. ++++++++++++++++++++++++++++++++YA LATIF++++++++++++++++++++++++++++++++++
.. +                                                                        +
.. + ScientiMate                                                            +
.. + Earth-Science Data Analysis Library                                    +
.. +                                                                        +
.. + Developed by: Arash Karimpour                                          +
.. + Contact     : www.arashkarimpour.com                                   +
.. + Developed/Updated (yyyy-mm-dd): 2017-07-01                             +
.. +                                                                        +
.. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scientimate.seqcolormap
=======================

.. code:: python

    CmapAKLines, CmapAKSet1, CmapAKPaired, CmapAKExtraLines, CmapAKTab10, CBrewerSet1, CBrewerPaired, Tab10 = scientimate.seqcolormap(dispout='no')

Description
-----------

Generate sequential colormap for drawing lines

Inputs
------

dispout='no'
    Define to display outputs or not ('yes': display, 'no': not display)

Outputs
-------

CmapAKLines
    7 Colors similar to matlab Lines in RGB (between 0 and 1)
CmapAKSet1
    9 Colors similar to colorbrewer2.org Set1 in RGB (between 0 and 1)
CmapAKPaired
    12 Colors similar to colorbrewer2.org Paired in RGB (between 0 and 1)
CmapAKExtraLines
    15 extra line colors similar to cbrewer qualitative Paired in RGB (between 0 and 1)
CmapAKTab10
    10 Colors similar to classic Tableau 10 color palettes in RGB (between 0 and 1)
CBrewerSet1
    9 Colors Set1 from Cynthia Brewer's ColorBrewer colorbrewer2.org in RGB (between 0 and 1)
CBrewerPaired
    12 Colors Paired from Cynthia Brewer's ColorBrewer colorbrewer2.org in RGB (between 0 and 1)
Tab10
    10 Colors from classic Tableau 10 color palettes in RGB (between 0 and 1)

Examples
--------

.. code:: python

    import scientimate as sm
    CmapAKLines,CmapAKSet1,CmapAKPaired,CmapAKExtraLines,CmapAKTab10,CBrewerSet1,CBrewerPaired,Tab10=sm.seqcolormap('yes')

References
----------

Colormap

* http://colorbrewer2.org
* http://matplotlib.org/cmocean/
* https://matplotlib.org/users/colormaps.html
* http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
* https://www.giss.nasa.gov/tools/panoply/colorbars/
* http://jdherman.github.io/colormap/

Source

* http://www.beautycolorcode.com  sample: http://www.beautycolorcode.com/0072bd  (RGB and hex code)
* http://www.htmlcsscolor.com/
* http://www.colorhexa.com/
* http://colorbrewer2.org/
* http://onlinehelp.tableau.com/current/pro/desktop/en-us/help.htm#formatting_create_custom_colors.html

.. License & Disclaimer
.. --------------------
..
.. Copyright (c) 2020 Arash Karimpour
..
.. http://www.arashkarimpour.com
..
.. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
.. IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
.. FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
.. AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
.. LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
.. OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
.. SOFTWARE.
