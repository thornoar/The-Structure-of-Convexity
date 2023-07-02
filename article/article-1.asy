if(!settings.multipleView) settings.batchView=false;
settings.tex="pdflatex";
settings.inlinetex=true;
deletepreamble();
defaultfilename="article-1";
if(settings.render < 0) settings.render=4;
settings.outformat="";
settings.inlineimage=true;
settings.embed=true;
settings.toolbar=false;
viewportmargin=(2,2);


usepackage("bm");
texpreamble("\def\V#1{\bm{#1}}");
defaultpen(fontsize(12pt));

import "LaTeX/Asymptote/General.asy" as general;



defaultpen(linewidth(.9pt));

picture part;

pair[] a = new pair[7];

a[0] = (-3, 0);
a[1] = (-1.5, 0);
a[2] = (-.3, 0);
a[3] = (1.5, .5);
a[4] = (3, .5);
a[5] = (1.5, -.5);
a[6] = (3, -.5);

pen dotpen = linewidth(4pt);
pen polytopedotpen = linewidth(5.1pt);
pen polytopepen = linewidth(1.1pt);

draw(part, a[6]--a[5], arrow = MidArcArrow(SimpleHead));
draw(part, a[4]--a[3], arrow = MidArcArrow(SimpleHead));
draw(part, a[1]--a[0], arrow = MidArcArrow(SimpleHead));

label(part, "\(a\)", a[1], 2*S);
label(part, "\(b\)", a[3], 2*N);
label(part, "\(c\)", a[5], 2*S);

void drawsegment (pair a, pair b, pen drawpen = currentpen)
{
draw(a--b, drawpen);
dot(a, dotpen);
dot(b, dotpen);
}

settings.render = 0;
// settings.outformat = "png";

size(9cm);

import three;
import graph3;

path3 myarc = rotate(18, Z) * Arc(c = O, normal = X, v1 = -Z, v2 = Z, n = 10);
surface backHemisphere = surface(myarc, angle1=0, angle2=180, c=O, axis = Z, n = 10);
surface frontHemisphere = surface(myarc, angle1 = 180, angle2 = 360, c = O, axis = Z, n = 10);
draw(backHemisphere, surfacepen = material(white+opacity(0.8), emissivepen=0.1*white), meshpen=gray(0.4));
draw(frontHemisphere, surfacepen=material(white+opacity(0.8),
emissivepen=0.1*white), meshpen=gray(0.4));



