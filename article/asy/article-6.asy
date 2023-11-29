if(!settings.multipleView) settings.batchView=false;
settings.tex="pdflatex";
defaultfilename="article-6";
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
size(10cm, 6cm);

pair a = (-2, 5), b = (-1, -3), x = (-4, 0), y = (4, 1.5), z = (0, 2);

label(a, "\(a\)", align = 1.1*N);
label(b, "\(b\)", align = S);
label(x, "\(x\)", align = W);
label(y, "\(y\)", align = E);
label(z, "\(z\)", align = .9*E+N);

drawsegment(a, y);
drawsegment(y, b);
drawsegment(x, b);
drawsegment(a, x);
drawsegment(a, b);
drawsegment(a, z);
drawsegment(z, b);
drawsegment(x, z);
drawsegment(z, y);
drawsegment(x, y);

label("\(1\)", a--x, align = (0,0), filltype = Fill(white));
label("\(1\)", a--y, align = (0,0), filltype = Fill(white));
label("\(1\)", x--b, align = (0,0), filltype = Fill(white));
label("\(1\)", y--b, align = (0,0), filltype = Fill(white));
label("\(2\)", a--z, align = (0,0), filltype = Fill(white));
label(Label("\(2\)", position = Relative(.25)), a--b, align = (0,0), filltype = Fill(white));
label("\(2\)", z--b, align = (0,0), filltype = Fill(white));
label(Label("\(2\)", position = Relative(.67)), x--y, align = (0,0), filltype = Fill(white));
label(Label("\(1\)", position = Relative(.35)), x--z, align = (0,0), filltype = Fill(white));
label("\(1\)", z--y, align = (0,0), filltype = Fill(white));
///>
