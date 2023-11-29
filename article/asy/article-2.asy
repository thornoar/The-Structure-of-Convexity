if(!settings.multipleView) settings.batchView=false;
settings.tex="pdflatex";
defaultfilename="article-2";
if(settings.render < 0) settings.render=4;
settings.outformat="";
settings.inlineimage=true;
settings.embed=true;
settings.toolbar=false;
viewportmargin=(2,2);

///<
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
///>

///<
size(7cm);

path convex = (-1,0)..(0,2)..(1,0)..(0,-1)..cycle;
path unconvex = (-1,0)..(0,1.5)..(1.5,.7)..(.2,0)..(.8,-.8)..(0,-1)..cycle;

draw(shift(-1,0)*convex);
draw(shift(2,0)*unconvex);

drawsegment((-1.25, 1.8), (-.11, .8));
drawsegment((-.9,-.8), (-1.66, 1.45));

label("convex", (-1, -1.6));

drawsegment((2.5, -.7), (3.1, .9), drawpen = currentpen+red);

label("not convex", (2, -1.6));
///>
