if(!settings.multipleView) settings.batchView=false;
settings.tex="pdflatex";
settings.inlinetex=true;
deletepreamble();
defaultfilename="article-3";
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

size(8cm);

draw((-.05,0)--(1,0), arrow = ArcArrow(SimpleHead));
draw((0,-.05)--(0,1), arrow = ArcArrow(SimpleHead));

pair v1 = (.7,.4);
pair v2 = (.7,1.1);
pair v3 = (1.8,1.1);
pair v4 = (1.8,.4);

real nudge = .1;
pair vnudge = (0, nudge);
pair hnudge = (nudge, 0);

draw(v1-vnudge--v1, dashed);
draw(v2+vnudge--v2, dashed);
draw(v2-hnudge--v2, dashed);
draw(v3+hnudge--v3, dashed);
draw(v3+vnudge--v3, dashed);
draw(v4-vnudge--v4, dashed);
draw(v4+hnudge--v4, dashed);
draw(v1-hnudge--v1, dashed);

filldraw(v1--v2--v3--v4--cycle, drawpen = blue+.8, fillpen = blue+opacity(.1));

dot(v1);
dot(v2);
dot(v3);
dot(v4, linewidth(4pt));

real rad = .02;

filldraw(circle(v1, rad), fillpen = red);
filldraw(circle(v3, rad), fillpen = red);

filldraw(circle(v2, rad), fillpen = green);
filldraw(circle((v1+v4)/2, rad), fillpen = green);
filldraw(circle((v3+v4)/2, rad), fillpen = green);
