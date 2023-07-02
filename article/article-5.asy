if(!settings.multipleView) settings.batchView=false;
settings.tex="pdflatex";
settings.inlinetex=true;
deletepreamble();
defaultfilename="article-5";
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



picture p, q, pq;

p.add(part);
q.add(part);
pq.add(part);

draw(p, a[5]{W}..{W}a[2], arrow = MidArcArrow(SimpleHead));
dot(p, a[1], polytopedotpen+red);
dot(p, a[3], polytopedotpen+red);
dot(p, a[5], dotpen);
draw(p, a[3]{W}..{W}a[2], p = polytopepen+red, arrow = MidArcArrow(SimpleHead));
draw(p, a[2]--a[1], p = polytopepen+red, arrow = MidArcArrow(SimpleHead));

label(p, "\(\dim P = 1\)", (0, -1.3));

add(shift(-6,0)*p);


draw(q, a[3]{W}..{W}a[2], arrow = MidArcArrow(SimpleHead));
dot(q, a[1], polytopedotpen+blue);
dot(q, a[5], polytopedotpen+blue);
dot(q, a[3], dotpen);
draw(q, a[5]{W}..{W}a[2], p = polytopepen+blue, arrow = MidArcArrow(SimpleHead));
draw(q, a[2]--a[1], p = polytopepen+blue, arrow = MidArcArrow(SimpleHead));

label(q, "\(\dim Q = 1\)", (0, -1.3));

add(q);


dot(pq, a[1], polytopedotpen+deepgreen);
dot(pq, a[3], polytopedotpen+deepgreen);
dot(pq, a[5], polytopedotpen+deepgreen);
draw(pq, a[5]{W}..{W}a[2], p = polytopepen+deepgreen, arrow = MidArcArrow(SimpleHead));
draw(pq, a[3]{W}..{W}a[2], p = polytopepen+deepgreen, arrow = MidArcArrow(SimpleHead));
draw(pq, a[2]--a[1], p = polytopepen+deepgreen, arrow = MidArcArrow(SimpleHead));

label(pq, "\(\dim \tp{P, Q} = 2\)", (0, -1.3));

add(shift(6,0)*pq);
size(483.69687pt,0,keepAspect=true);
