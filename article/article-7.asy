if(!settings.multipleView) settings.batchView=false;
settings.tex="pdflatex";
settings.inlinetex=true;
deletepreamble();
defaultfilename="article-7";
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



pair[] b = new pair[6];

for (int i = 0; i < b.length; ++i) {b[i] = (i, 0);}

for (int i = b.length-1; i > 1; --i)
{
draw(b[i] -- b[i-1], arrow = MidArcArrow(SimpleHead));
}

dot(b[0], dotpen);
label("\(y\)", b[0], align = 1.5*S);

for (int i = 1; i < b.length; ++i)
{
dot(b[i], dotpen);
label((string)i, b[i], align = 1.5*S);
}

label("\ldots", (b.length-.6, 0));

dot("\(\infty_1\)", align = 1.5*E, (b.length+.5, .5));
dot("\(\infty_2\)", align = 1.5*E, (b.length+.5, -.5));

draw((b.length+.5, .5) -- (b.length-.5, 0), arrow = MidArcArrow(SimpleHead), margin = Margin(begin = 0, end = 5));
draw((b.length+.5, -.5) -- (b.length-.5, 0), arrow = MidArcArrow(SimpleHead), margin = Margin(begin = 0, end = 5));

// drawWFunction(l = new string[]{"1", "2"}, mode = -1);
size(483.69687pt,0,keepAspect=true);
