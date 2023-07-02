settings.outformat = "pdf";
settings.render = 0;

size(841mm, 1189mm);

// import "LaTeX/General.asy" as general;
import "smoothmanifold/smoothmanifold.asy" as smooth;

real w = 300;
real h = w * 1189 / 841; //h = 424.13793103448273, h/2 = 212
real w2 = w/2;
real h2 = h/2;
real margin = 7;

draw((-w2, -h2)--(w2, h2), invisible);
draw((-w2, h2)--(w2, -h2), invisible);

// ? Margins and control points

pair wn = (-w2 + margin, h2 - margin);
pair en = (w2 - margin, h2 - margin);
pair es = (w2 - margin, -h2 + margin);
pair ws = (-w2 + margin, -h2 + margin);

path outerframe = (-w2, h2)--(w2, h2)--(w2, -h2)--(-w2, -h2)--cycle;

path innerframe = wn--en--es--ws--cycle;

real trailw = 7;
real blockrw = (w - 2*margin)/3;

real bodyu = 160;
real bodymid = -135;
real midhpoint = 20;
real righthpoint = -80;
real bodyd = -188;

real c1 = -w2+margin;
real c2 = c1+blockrw;
real c3 = c2+blockrw;
real c4 = w2-margin;

real convexl = c1;
real convexr = c2 - trailw/2;
real convexu = bodyu;
real convexd = bodymid + trailw/2;
pair convex1 = (convexl, convexu);
pair convex2 = (convexr, convexu);
pair convex3 = (convexr, convexd);
pair convex4 = (convexl, convexd);
pair convex12 = (convex1 + convex2)/2;
pair convex23 = (convex2 + convex3)/2;
pair convex34 = (convex3 + convex4)/2;
pair convex41 = (convex4 + convex1)/2;

real orderl = c2 + trailw/2;
real orderr = c3 - trailw/2;
real orderu = bodyu;
real orderd = midhpoint + trailw/2;
pair order1 = (orderl, orderu);
pair order2 = (orderr, orderu);
pair order3 = (orderr, orderd);
pair order4 = (orderl, orderd);
pair order12 = (order1 + order2)/2;
pair order23 = (order2 + order3)/2;
pair order34 = (order3 + order4)/2;
pair order41 = (order4 + order1)/2;

real metricl = orderl;
real metricr = orderr;
real metricu = midhpoint - trailw/2;
real metricd = bodymid + trailw/2;
pair metric1 = (metricl, metricu);
pair metric2 = (metricr, metricu);
pair metric3 = (metricr, metricd);
pair metric4 = (metricl, metricd);
pair metric12 = (metric1 + metric2)/2;
pair metric23 = (metric2 + metric3)/2;
pair metric34 = (metric3 + metric4)/2;
pair metric41 = (metric4 + metric1)/2;

real manl = c3 + trailw/2;
real manr = c4;
real manu = bodyu;
real mand = righthpoint + trailw/2;
pair man1 = (manl, manu);
pair man2 = (manr, manu);
pair man3 = (manr, mand);
pair man4 = (manl, mand);
pair man12 = (man1 + man2)/2;
pair man23 = (man2 + man3)/2;
pair man34 = (man3 + man4)/2;
pair man41 = (man4 + man1)/2;

real refl = manl;
real refr = manr;
real refu = righthpoint - trailw/2;
real refd = bodymid + trailw/2;
pair ref1 = (refl, refu);
pair ref2 = (refr, refu);
pair ref3 = (refr, refd);
pair ref4 = (refl, refd);
pair ref12 = (ref1 + ref2)/2;
pair ref23 = (ref2 + ref3)/2;
pair ref34 = (ref3 + ref4)/2;
pair ref41 = (ref4 + ref1)/2;

real addl = c1;
real addr = c4;
real addu = bodymid - trailw/2;
real addd = bodyd;
pair add1 = (addl, addu);
pair add2 = (addr, addu);
pair add3 = (addr, addd);
pair add4 = (addl, addd);
pair add12 = (add1 + add2)/2;
pair add23 = (add2 + add3)/2;
pair add34 = (add3 + add4)/2;
pair add41 = (add4 + add1)/2;

// ? bubbles

real inonein = 9.09091;
real boxtitleinches = .5;
real textinches = .3;

real boxtitleheight = inonein * boxtitleinches;
real textheight = inonein * textinches;
real titleboxoffset = 4;
real textoffset = 2;

pen bubblepen = linewidth(2mm);
pen boxtitlebubblepen = linewidth(1.5mm);
pen textbubblepen = linewidth(1mm);
pen arrowpen = linewidth(.7mm);
pen dashedpen = linewidth(.5mm)+dashed;
pen dotpen = linewidth(2.5mm);
pen p = linewidth(.7mm);

path boxtitlebubble (string s, pair a, real offset = titleboxoffset)
{
    picture sub;

    label(sub, s);

    real r = size(sub).x/size(sub).y;

    if(r > 15) r+=4;

    real th = boxtitleheight;
    real tw = r*th;

    real offsetrel = .7;

    pair p1 = a+(0, -th/2-offset*(offsetrel+.1));
    pair p2 = a+(-tw/2+offset/2, -th/2-offset*offsetrel);
    pair p3 = a+(-tw/2-offset*1.4, 0);
    pair p4 = a+(-tw/2+offset/2, th/2+offset*offsetrel*1.08);
    pair p5 = a+(0, th/2+offset*(offsetrel+.1));
    pair p6 = a+(tw/2-offset/2, th/2+offset*offsetrel);
    pair p7 = a+(tw/2+offset*1.3, 0);
    pair p8 = a+(tw/2-offset/2, -th/2-offset*offsetrel);

    // dot(p1, red+10);
    // dot(p2, red+10);
    // dot(p3, red+10);
    // dot(p4, red+10);
    // dot(p5, red+10);
    // dot(p6, red+10);
    // dot(p7, red+10);
    // dot(p8, red+10);

    return p1{W}..p2{W+.1*N}..p3{N+.2*E}..p4{E+.1*N}..p5{E}..p6{E+.1*S}..p7{S+.2*E}..p8{W+.1*S}..cycle;
}
path textbubble (string s, pair a, real offset = textoffset)
{
    picture sub;

    label(sub, s);

    real r = size(sub).x/size(sub).y;

    real udot = .02;
    real udot2 = .01;
    real offsetrel = .8;

    if(r > 10) 
    {
        r+=6.5;
        udot /= 1.1;
        // udot2 /= 8;
        offsetrel = 1.1;
    }

    real th = textheight;
    real tw = r*th;
    real corneroffsetm = tw*udot2;

    pair p1 = a+(0, -th/2-offset*(offsetrel+tw*udot));
    pair p2 = a+(-tw/2+offset*corneroffsetm, -th/2-offset*(offsetrel+.02));
    pair p3 = a+(-tw/2-offset*1.4, 0);
    pair p4 = a+(-tw/2+offset*corneroffsetm, th/2+offset*offsetrel*1.08);
    pair p5 = a+(0, th/2+offset*(offsetrel+tw*udot));
    pair p6 = a+(tw/2-offset*corneroffsetm, th/2+offset*offsetrel);
    pair p7 = a+(tw/2+offset*1.3, 0);
    pair p8 = a+(tw/2-offset*corneroffsetm, -th/2-offset*offsetrel);

    return p1{W}..p2{W+.2*N}..p3{N+.2*W}..p4{E+.2*N}..p5{E}..p6{E+.2*S}..p7{S+.2*E}..p8{W+.3*S}..cycle;
}
void addboxtitle (string s, pair a, pen contourpen = boxtitlebubblepen, real inches = boxtitleinches, real offset = titleboxoffset)
{
    string label = "\resizebox{!}{" + (string)inches + "in}{\scshape " + s + "}";

    filldraw(boxtitlebubble(s = label, a = a, offset = offset), fillpen = white, drawpen = contourpen);
    label(label, a);
}
void addtext (string s, pair a, pen contourpen = textbubblepen, real inches = textinches, real offset = textoffset)
{
    string label = "\resizebox{!}{" + (string)inches + "in}{\scshape " + s + "}";

    filldraw(textbubble(s = label, a = a, offset = offset), fillpen = white, drawpen = contourpen);
    label(label, a);
}
void plainaddtext (string s, pair a, pen contourpen = textbubblepen, real inches = textinches, real offset = textoffset)
{
    string label = "\resizebox{!}{" + (string)inches + "in}{" + s + "}";

    filldraw(textbubble(s = label, a = a, offset = offset), fillpen = white, drawpen = contourpen);
    label(label, a);
}
void drawsegment (picture pic = currentpicture, pair a, pair b, pen drawpen = textbubblepen, pen dotpen = dotpen)
{
    draw(pic, a--b, drawpen);
    dot(pic, a, dotpen);
    dot(pic, b, dotpen);
}
void drawpolytope (picture pic = currentpicture, pair[] v, pen drawpen = textbubblepen, pen dotpen = dotpen)
{
    for (int i = 0; i < v.length; ++i)
    {
        draw(pic, v[i]--v[(i+1)%v.length], drawpen);
        dot(pic, v[i], dotpen);
    }
}
path arraypath (pair[] v)
{
    path res = v[0]--v[1];
    for (int i = 2; i < v.length; ++i)
    {
        res = res -- v[i];
    }
    res = res--cycle;
    return res;
}
void filldrawpolytope (picture pic = currentpicture, pair[] v, pen drawpen = textbubblepen, pen dotpen = dotpen, pen fillpen = cyan)
{
    fill(pic, arraypath(v), p = fillpen);
    drawpolytope(pic, v, drawpen = drawpen, dotpen = dotpen);
}

path titlebubble = shift(0,3)*((0, 166){W}..(-50, 166)..(-114, 169){W+.13*N}..(-125, 178.5){N+.05*E}..(-112.5, 194){E+.15*N}..{E}(0, 197)..{E+.15*S}(115, 193.5)..(126, 184.5){S+.05*W}..(114, 169.5){W+.1*S}..cycle);

path convexbubble = (
    (convex12+(0,1)){E}..(convex2+(-13,0)){E+.11*S}..{S+.01*E}(convex2+(0,-25))..
    convex23{S}..(convex3+(0,17)){S+.01*W}..{W+.05*S}(convex3+(-16,0))..
    (convex34+(0,-.5)){W}..(convex4+(17,0)){W+.08*N}..{N+.01*W}(convex4+(0,15))..
    convex41{N}..(convex1+(0,-30)){N}..{E+.1*N}(convex1+(18,0))..cycle
);


path orderbubble = (
    (order12+(0,.6)){E}..(order2+(-20,0)){E+.07*S}..{S+.05*E}(order2+(0,-20))..
    order23{S}..(order3+(0,17)){S+.02*W}..{W+.07*S}(order3+(-16,0))..
    (order34+(0,-.5)){W}..(order4+(16,0)){W+.1*N}..{N+.03*W}(order4+(0,13))..
    order41{N}..(order1+(0,-19)){N+.04*E}..{E+.1*N}(order1+(17,0))..cycle
);


path metricbubble = (
    (metric12+(0,.7)){E}..(metric2+(-16,0)){E+.12*S}..{S+.02*E}(metric2+(0,-18))..
    metric23{S}..(metric3+(0,15)){S+.03*W}..{W+.07*S}(metric3+(-16,0))..
    (metric34+(0,-.6)){W}..(metric4+(17,0)){W+.07*N}..{N+.04*W}(metric4+(0,15))..
    metric41{N}..(metric1+(0,-16)){N+.03*E}..{E+.1*N}(metric1+(18,0))..cycle
);


path manbubble = (
    (man12+(0,1)){E}..(man2+(-23,0)){E+.11*S}..{S+.05*E}(man2+(0,-28))..
    man23{S}..(man3+(0,17)){S+.01*W}..{W+.05*S}(man3+(-16,0))..
    (man34+(0,-.5)){W}..(man4+(17,0)){W+.08*N}..{N+.01*W}(man4+(0,15))..
    man41{N}..(man1+(0,-30)){N}..{E+.1*N}(man1+(18,0))..cycle
);


path refbubble = (
    (ref12+(0,.6)){E}..(ref2+(-20,0)){E+.07*S}..{S+.05*E}(ref2+(0,-20))..
    ref23{S}..(ref3+(0,17)){S+.02*W}..{W+.07*S}(ref3+(-16,0))..
    (ref34+(0,-.5)){W}..(ref4+(16,0)){W+.1*N}..{N+.03*W}(ref4+(0,13))..
    ref41{N}..(ref1+(0,-19)){N+.04*E}..{E+.1*N}(ref1+(17,0))..cycle
);


path addbubble = (
    (add12+(0,.7)){E}..(add2+(-16,0)){E+.05*S}..{S+.02*E}(add2+(0,-18))..
    add23{S}..(add3+(0,24)){S+.05*W}..{W+.07*S}(add3+(-24,0))..
    (add34+(0,-1.8)){W}..(add4+(20,0)){W+.07*N}..{N+.04*W}(add4+(0,25))..
    add41{N}..(add1+(0,-17)){N+.03*E}..{E+.04*N}(add1+(17,0))..cycle
);

// ? the drawing

fill(outerframe ^^ reverse(titlebubble) ^^ reverse(convexbubble) ^^ reverse(orderbubble) ^^ reverse(metricbubble) ^^ reverse(manbubble) ^^ reverse(refbubble) ^^ reverse(addbubble), cyan+opacity(.8));

real opacity = .92;

fill(convexbubble, springgreen+lightcyan+opacity(opacity));
fill(orderbubble, magenta+cyan+opacity(opacity));
fill(metricbubble, purple+cyan+opacity(opacity));
fill(manbubble, mediumblue+cyan+opacity(opacity));
fill(refbubble, deepcyan+opacity(opacity));
fill(addbubble, yellow+cyan+opacity(opacity));

draw(titlebubble, bubblepen);
draw(convexbubble, bubblepen);
draw(orderbubble, bubblepen);
draw(metricbubble, bubblepen);
draw(manbubble, bubblepen);
draw(refbubble, bubblepen);
draw(addbubble, bubblepen);

label((0, 185), "{\resizebox{25in}{!}{\scshape The Structure of Convexity}}");

addboxtitle(s = "Internal Theory", a = ((convexl+convexr)/2+1.5, 151));
addboxtitle(s = "Inducing Structure", a = ((orderl+orderr)/2, 151));
addboxtitle(s = "Induced Structure", a = ((manl+manr)/2, 150));
addboxtitle(s = "References", a = ((refl+refr)/2, -91.2));
addboxtitle(s = "Uniquely Geodesic Spaces", a = (-85, -147));

// ? filling the internal theory

path itcont1 = shift(-97, 120)*rotate(165)*scale(7.7)*samplesmooth(2).contour;
filldraw(itcont1, fillpen = white, drawpen = boxtitlebubblepen);

addtext(s = "Idea", a = (-118, 94));
addtext(s = "Definition", a = (-94, 135));

picture sub1;

path convex = (-1,0)..(0,2)..(1,0)..(0,-1)..cycle;
path unconvex = shift(-.4, 0)*((-1,0)..(0,1.5)..(1.5,.7)..(.2,0)..(.8,-.8)..(0,-1)..cycle);

draw(sub1, shift(-1,0)*convex, textbubblepen);
draw(sub1, shift(2,0)*unconvex, textbubblepen);

drawsegment(sub1, (-1.25, 1.8), (-.11, .8), drawpen = textbubblepen);
drawsegment(sub1, (-.9,-.8), (-1.66, 1.45), drawpen = textbubblepen);

drawsegment(sub1, (2.1, -.7), (2.7, .9), drawpen = textbubblepen+red);

add(shift(-118, 108.5)*scale(6.5)*sub1);

label("\Huge \bfseries \((X, \mathcal{C})\) --- \textit{convex space:}", (-88, 127));
label("\Huge \bfseries (1) \(\emptyset, \ X \in \mathcal{C}\)", (-94, 122), align = E);
label("\Huge \bfseries (2) \(\mathcal{A} \subset \mathcal{C} \Rightarrow \cap \mathcal{A} \in \mathcal{C}\)", (-94, 117), align = E);
label("\Huge \bfseries (3) \(\mathcal{N} \subset \mathcal{C} \Rightarrow \cup \mathcal{N} \in \mathcal{C}\)", (-94, 112), align = E);


path itcont2 = shift(-75, 70)*rotate(-55)*scale(11)*samplesmooth(1).contour;
filldraw(itcont2, fillpen = white, drawpen = boxtitlebubblepen);

addtext(s = "Polytope", a = (-81, 94));

drawsegment((-94, 75), (-87, 85));
pair q1 = (-79, 84);
pair q2 = (-70, 88);
pair q3 = (-63, 79);
pair q4 = (-67, 73);
// fill(q1--q2--q3--q4--cycle, mediumblue);
filldrawpolytope(new pair[]{q1,q2,q3,q4}, fillpen = mediumblue);

pair tet1 = (-90, 69);
pair tet2 = (-82, 79);
pair tet3 = (-74, 70);
pair tet4 = (-80, 65.5);
fill(tet1--tet2--tet3--tet4--cycle, deepgreen);
fill(tet1--tet2--tet4--cycle, heavygreen);
drawpolytope(new pair[]{tet1, tet2, tet3, tet4});
draw(tet1--tet3, dashedpen);
draw(tet2--tet4, textbubblepen);

addtext(s = "Dimension", a = (-77, 58));

picture sub2;
int n = 7;
pair[] reg = new pair[n];
for (int i = 0; i < n; ++i)
{
    reg[i] = (cos(2*pi*i/n)-.2, sin(2*pi*i/n));
}
// path regp = arraypath(reg);
// fill(sub2, regp, mediumyellow);
filldrawpolytope(sub2, reg, fillpen = mediumyellow);
// label(sub2, "\Huge \(P\)", reg[floor(n/2)], align = 2*W+S);
pair tt1 = (-3, 1.5);
pair tt2 = (2.6, 1.05);
pair tt3 = (-.5, -2.1);
drawpolytope(sub2, new pair[]{tt1,tt2,tt3}, drawpen = dashedpen);
add(shift(-76, 42)*scale(5.7)*sub2);

path itcont3 = shift(-95, 25)*scale(16)*rotate(-70)*samplesmooth(2).holes[0].contour;
filldraw(itcont3, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Freedom", a = (-119, 37.5));
pair f1 = (-130, 61);
pair f2 = (-121, 73);
pair f3 = (-111, 63);
pair f4 = (-113, 51);
pair f5 = (-124, 47);
pair[] f = new pair[]{f1,f2,f3,f4,f5};
// fill(f1--f2--f3--f4--f5--cycle, royalblue+opacity(.4));
filldrawpolytope(f, drawpen = dashedpen, fillpen = royalblue+opacity(.4));
for (int i = 0; i < f.length; ++i)
{
    draw(circle(f[i], 1.5), 2+red);
}

path itcont4 = shift(-90, -7)*scale(7.6)*rotate(170)*((-4.5,1)..(-2.3,8.5)..(1.7,1.8)..(5.6,1)..(2,-4)..(-3,-2.5)..cycle);
filldraw(itcont4, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Hyperplane", a = (-102, 19.7));
label("\Huge \textit{(Maximal net of polytopes of same dim.)}", (-100, 12.5));
pair coordc = (-117, -1);
real xl = 21;
real yl = 19;
real zl = 13;
real zr = .7;
real y = zl/(2*sqrt(zr^2 + 1));
real x = zr*y;
draw((coordc+(-xl/2, 0))--(coordc+(xl/2, 0)), arrow = ArcArrow(SimpleHead), dashedpen);
draw((coordc+(0, -yl/2))--(coordc+(0, yl/2)), arrow = ArcArrow(SimpleHead), dashedpen);
draw((coordc+(-x, -y))--(coordc+(x, y)), arrow = ArcArrow(SimpleHead), dashedpen);
pair par1 = (-110, -9);
pair par2 = (-75, 8);
pair par3 = (-62, -7);
pair par4 = par1+(par3-par2);
fill(par1--par2--par3--par4--cycle, grey+opacity(.8));
filldrawpolytope(new pair[]{(-91, -17), (-100, -14), (-96, -5), (-80, 2), (-73, -2), (-71, -8)}, fillpen = lightmagenta+mediumblue);
filldrawpolytope(new pair[]{(-96, -10), (-88, -7), (-91, -14)}, fillpen = purple+mediumcyan);
filldrawpolytope(new pair[]{(-83, -10), (-84, -2.5), (-78.5, -1), (-75.5, -8)}, fillpen = blue+mediumgreen);
pair arrowpos1 = (-100, -9.5);
draw(arrowpos1--(arrowpos1+5*dir(175)), p = arrowpen, arrow = ArcArrow(SimpleHead));
pair arrowpos2 = (-77, 2);
draw(arrowpos2--(arrowpos2+4.1*dir(60)), p = arrowpen, arrow = ArcArrow(SimpleHead));
pair arrowpos3 = (-70, -6);
draw(arrowpos3--(arrowpos3+5*dir(-7)), p = arrowpen, arrow = ArcArrow(SimpleHead));
pair arrowpos4 = (-94.5, -18);
draw(arrowpos4--(arrowpos4+4.1*dir(-115)), p = arrowpen, arrow = ArcArrow(SimpleHead));
addtext(s = "TPUL", a = (-78, -68));
pair[] u1 = new pair[]{(-84, -58.5), (-90, -47), (-76, -36), (-65, -48), (-70, -60)}; 
pair[] u2 = new pair[]{(-77, -53), (-91, -38), (-71, -24), (-63, -36)}; 
pair in1 = intersectionpoints(arraypath(u1), arraypath(u2))[0];
pair in2 = intersectionpoints(arraypath(u1), arraypath(u2))[1];
fill(u1[0]--u1[1]--in1--u2[0]--in2--u1[3]--u1[4]--cycle, p = mediumblue);
fill(u2[1]--u2[2]--u2[3]--in2--u1[2]--in1--cycle, p = mediumyellow);
drawpolytope(u1);
drawpolytope(u2);
filldrawpolytope(new pair[]{(-79, -48.5), (-81, -42), (-72, -44)}, fillpen = mediumgreen);
draw(u1[1]--u2[1], dashedpen);
draw(u1[3]--u2[3], dashedpen);

path itcont5 = shift(-121, -51)*yscale(.9)*xscale(1.05)*scale(9)*rotate(122)*samplesmooth(1).contour;
filldraw(itcont5, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Subsets", a = (-117.5, -26));
path subset = shift(-118, -48)*scale(3.2)*rotate(0)*xscale(.8)*((0,5){E}..(3.5,2)..(3,-4)..(0,-3){W}..(-3,-4)..(-3.5,2)..cycle);
draw(subset, textbubblepen);
pair x1 = (-117, -42);
dot(x1, dotpen);
label("\Huge \(x\)", x1, align = 2.5*W+S);
drawpolytope(new pair[]{(-123, -45), (-120, -39), (-113, -40), (-112, -45)}, drawpen = dashedpen, dotpen = linewidth(2mm));
label("\Huge \(\eta_x\)", (-123, -45), align = 2.2*S+2*E);
pair y1 = point(subset, 3.5);
dot(y1, dotpen);
label("\Huge \(y\)", y1, align = 2*E+1.4*S);
drawpolytope(new pair[]{y1, y1+(-4, 4), y1+(3, 7)}, drawpen = dashedpen, dotpen = linewidth(2mm));
label("\Huge \(\eta_y\)", y1+(3,7), align = 2.5*S+1.7*E);
label("\Huge \(\eta_x = \eta_y = \mathrm{dim}\ A\)", (-116, -66));

pair[] shift (pair[] v, pair sh)
{
    return sequence(new pair (int i){return v[i]+sh;}, v.length);
}

path itcont6 = shift(-102.5, -98)*scale(11.5)*rotate(30)*samplesmooth(1).contour;
filldraw(itcont6, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Theorem", a = (-100, -84));
label("\Huge \textit{TPUL is equivalent to:}", (-100, -92));
picture sub3;
filldrawpolytope(sub3, new pair[]{(-91, -17), (-100, -14), (-96, -5), (-80, 2), (-73, -2), (-71, -8)}, fillpen = lightmagenta+mediumblue);
filldrawpolytope(sub3, new pair[]{(-96, -10), (-88, -7), (-91, -14)}, fillpen = purple+mediumcyan);
filldrawpolytope(sub3, new pair[]{(-83, -10), (-84, -2.5), (-78.5, -1), (-75.5, -8)}, fillpen = blue+mediumgreen);
add(shift(-118, -105)*rotate(-90)*xscale(1.5)*rotate(70)*scale(.9)*shift((88, 7))*sub3);
picture sub4;
pair[] pentagon = new pair[]{(-2, -2), (-3, 1), (1, 3), (3, 0), (1, -3)};
drawpolytope(sub4, pentagon);
picture sub5;
fill(sub5, arraypath(rotate(150)*scale(3.9)*pentagon), magenta+lightcyan);
fill(sub5, arraypath(rotate(80)*scale(2.7)*pentagon), mediumgreen+mediumcyan);
fill(sub5, arraypath(rotate(20)*scale(1.7)*pentagon), mediumblue+mediumcyan);
fill(sub5, arraypath(pentagon), lightyellow);
sub5.add(sub4);
sub5.add(rotate(20)*scale(1.7)*sub4);
sub5.add(rotate(80)*scale(2.7)*sub4);
sub5.add(rotate(150)*scale(3.9)*sub4);
add(shift(-79, -106)*sub5);
label("\resizebox{1in}{!}{\(\Longleftrightarrow\)}", (-96.5, -105));


// ? filling the inducing structures

void drawchain (picture pic = currentpicture, pair[] v, pen drawpen = linewidth(.7mm), pen dotpen = dotpen, arrowbar arrow = MidArcArrow(SimpleHead))
{
    for (int i = 0; i < v.length; ++i)
    {
        dot(pic, v[i], dotpen);
        if (i < v.length-1) {draw(pic, v[i+1]--v[i], p = drawpen, arrow = arrow);}
    }
}

path iscont1 = shift(-1, 119)*yscale(.89)*scale(1.4)*rotate(-90)*shift(-path_middle(itcont5))*itcont5;
filldraw(iscont1, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Order", a = (-3, 138));
picture sub6;
drawchain(sub6, new pair[]{(0,0), (2,-1), (4,-1), (6,-1), (8,0), (10,0)});
drawchain(sub6, new pair[]{(0,-2), (2,-1)});
drawchain(sub6, new pair[]{(6,-1), (8,-2), (10,-2), (12,-2), (14,-1), (16,-1)});
drawchain(sub6, new pair[]{(12, -2), (14,-3)});
drawchain(sub6, new pair[]{(-1,-4), (1,-4), (3,-3), (5,-3)});
drawchain(sub6, new pair[]{(9,-4), (11,-5), (13,-5)});
drawchain(sub6, new pair[]{(5,-3), (7,-4), (9,-4), (11,-3)}, drawpen = red+linewidth(.8mm), dotpen = red+linewidth(3mm));
drawchain(sub6, new pair[]{(3,-5), (5,-5), (7,-4)}, drawpen = red+linewidth(.8mm), dotpen = red+linewidth(3mm));
drawchain(sub6, new pair[]{(0,-6), (2,-6), (4,-7), (6,-7), (8,-6), (10,-6), (12,-7)});
real rad = 1.5/3.5;
draw(sub6, circle((5,-3), rad), 2+blue);
draw(sub6, circle((3,-5), rad), 2+blue);
draw(sub6, circle((11,-3), rad), 2+blue);
add(shift(-26, 126)*scale(3.5)*sub6);
label("\Huge \textbf{Theorem:} \textit{all order convexities are free.}", (0, 131));

path iscont2 = shift(-2, 72)*rotate(-10)*xscale(.9)*yscale(.9)*reflect((0,0), (0,1))*shift(-path_middle(itcont1))*shift(-97, 120)*rotate(165)*scale(7.7)*((-5,-1.2)..(-.6,.8)..(1.5, 6.2)..(0.5,-2.5)..(-1.5,-2.7)..cycle);
filldraw(iscont2, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Linearity", a = (-15, 89));
picture sub7;
draw(sub7, (0, 0)--(10,0), textbubblepen);
dot(sub7, Label("\Huge \(a\)", align = 3*S), (2,0), dotpen);
dot(sub7, Label("\Huge \(x\)", align = 3*S), (4.5,0), dotpen);
dot(sub7, Label("\Huge \(y\)", align = 3*S), (7,0), dotpen);
dot(sub7, Label("\Huge \(b\)", align = 3*S), (9,0), dotpen);
draw(sub7, (7,0)--(4.5,0), textbubblepen, arrow = MidArcArrow(SimpleHead));
add(shift(-30, 80)*scale(5)*sub7);
addtext(s = "1-Affinity", a = (18.5, 42));
picture sub8;
path set = rotate(80)*scale(5.6)*yscale(.7)*rotate(-50)*samplesmooth(1).contour;
draw(sub8, set, textbubblepen);
pair x2 = (-4,-5);
pair y2 = (5,1);
draw(sub8, (x2+.9*(x2-y2))--(y2+1*(y2-x2)), textbubblepen);
dot(sub8, Label("\Huge \(x\)", align = 2*E+2*S),x2, dotpen);
dot(sub8, Label("\Huge \(y\)", align = 2*E+2*S), y2, dotpen);
add(shift(17.5, 63)*sub8);

path iscont3 = shift(-20, 49)*yscale(.75)*shift(-path_middle(itcont5))*itcont5;
filldraw(iscont3, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "\(n\)-Affinity", a = (-20, 62));
picture sub9;
draw(sub9, (-4,3)--(4,3)--(4,-3)--(-4,-3)--cycle, textbubblepen);
draw(sub9, (-4, -2)--(4,0), textbubblepen);
label(sub9, "\Huge \(L\)", (-2, -1.5), align = 2*S+E);
filldrawpolytope(sub9, new pair[]{(-3,0), (-2.3, 2), (-.4, .3)}, fillpen = cyan+lightblue);
label(sub9, "\Huge \(A\)", (4,3), align = 4*(W+S));
label(sub9, "\Huge \(B\)", (4,-3), align = 4*(W+N));
drawsegment(sub9, (1.6, .8), (.6, -2.1), drawpen = red+textbubblepen);
draw(sub9, circle(intersectionpoints((-4, -2)--(4,0), (1.6, .8)--(.6, -2.1))[0], 1.5/3.3), blue+2);
add(shift(-20, 46)*scale(3.3)*sub9);

path iscont4 = shift(1, -11)*yscale(.95)*rotate(73)*scale(1.5)*shift(-path_middle(itcont3))*itcont3;
filldraw(iscont4, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Metric", a = (-3, 8));
picture sub10;
drawsegment(sub10, (0,0), (8,-2));
label(sub10, "\Huge \(a\)", (0,0), align = 2*W+S);
label(sub10, "\Huge \(b\)", (8,-2), align = W+2*S);
dot(sub10, Label("\Huge \(x\)", align = 2*W+2*S), .7*(8,-2), dotpen);
label(sub10, "\Huge \(d(a, b) = d(a, x) + d(x, b)\)", (3.5, -4.25));
add(shift(-26, 2.5)*scale(3)*sub10);
picture sub11;
path convex = (-1,0)..(0,2)..(1,0)..(0,-1)..cycle;
filldraw(sub11, shift(-1,0)*convex, drawpen = textbubblepen, fillpen = mediumyellow);
drawsegment(sub11, (-1.25, 1.8), (-.11, .8), drawpen = textbubblepen);
drawsegment(sub11, (-.9,-.8), (-1.66, 1.45), drawpen = textbubblepen);
add(shift(32,-20)*scale(11)*rotate(0)*sub11);
addtext(s = "Segmential", a = (-7, -24));

path iscont5 = shift(0,-50)*yscale(.7)*rotate(145)*shift(-path_middle(itcont1))*itcont1;
filldraw(iscont5, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Finite", a = (10, -64.5));
pair fq1 = (-3,1);
pair fq2 = (1,3);
pair fq3 = (3,-1);
pair fq4 = (-2,-3);
picture sub12;
dot(sub12, fq1, dotpen);
dot(sub12, fq2, dotpen);
dot(sub12, fq3, dotpen);
dot(sub12, fq4, dotpen);
picture sub13;
drawpolytope(sub13, new pair[]{fq1,fq2,fq3,fq4});
draw(sub13, fq1--fq3, p);
draw(sub13, fq2--fq4, p);
picture sub14;
fill(sub14, fq1--fq2--fq3--fq4--cycle, lightcyan+magenta);
sub14.add(sub13);
draw(sub14, (.2*fq1+.8*fq2)--(.3*fq3+.7*fq4), p = p);
draw(sub14, (.85*fq1+.15*fq2)--(.4*fq2+.6*fq3), p = p);
draw(sub14, (.25*fq1+.75*fq4)--(.65*fq2+.35*fq3), p = p);
pair pos1 = (-25, -43);
pair pos2 = (0, -52);
pair pos3 = (27, -51);
real sc = 2.5;
add(shift(pos1)*scale(sc)*sub12);
add(shift(pos2)*scale(sc)*sub13);
add(shift(pos3)*scale(sc)*sub14);
real mg = 22;
real curve = .2;
draw(curved_path(pos1, pos2, curve = curve), arrow = ArcArrow(SimpleHead), margin = Margin(mg), p = p);
draw(curved_path(pos2, pos3, curve = curve), arrow = ArcArrow(SimpleHead), margin = Margin(mg), p = p);

path iscont6 = shift(-21,-96)*scale(1)*rotate(3)*reflect((0,0), (0,1))*shift(-path_middle(itcont5))*itcont5;
filldraw(iscont6, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Join", a = (-21, -73));
picture sub15;
pair x3 = (0,0);
pair j1 = (4,3);
pair j2 = (6,2);
pair j3 = (7,-1);
pair j4 = (4,-3);
path jp = j1--j2--j3--j4;
int n = 30;
fill(sub15, arraypath(new pair[]{j1,j2,j3,j4}), green+cyan);
for (int i = 0; i < n; ++i)
{
    draw(sub15, x3--point(jp, i/n*length(jp)), darkgrey+linewidth(.3mm));
}
draw(sub15, x3--j1, p);
draw(sub15, x3--j2, p);
draw(sub15, x3--j3, p);
draw(sub15, x3--j4, p);
dot(sub15, Label("\Huge \(x\)", align = 2*S+W), x3, dotpen);
drawpolytope(sub15, new pair[]{j1,j2,j3,j4});
label(sub15, "\Huge \(S\)", (j3+j4)/2, align = 2*S+E);
label(sub15, "\Huge \(x \vee S\)", (x3+j4)/2, align = 3*W+5*S);
add(shift(-34, -91)*scale(4)*sub15);
label("\Huge \textit{Join-commutative:}", (-21, -108.5));
label("\Huge \(x \vee \langle F \rangle = \langle x \cup F \rangle\)", (-21, -114));

path iscont7 = shift(20.5, -100)*reflect((0,0),(0,1))*xscale(.95)*yscale(1.28)*shift(-path_middle(iscont3))*iscont3;
filldraw(iscont7, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Theorem", a = (21, -82));

real ystep = -5;
real xpos = 21;
real ypos = -92;
label("\Huge \textit{Finite-segmential}", (xpos, ypos));
label("\Huge \textit{2-Affine}", (xpos, ypos+ystep));
label("\Huge \textit{TPUL}", (xpos, ypos+2*ystep));
label("\resizebox{!}{.4in}{\(\Downarrow\)}", (xpos, ypos+3.5*ystep));
label("\Huge \textit{Free}", (xpos, ypos+5*ystep));


// ? filling the uniquely geodesic spaces

path ugscont1 = shift((-81,-171))*xscale(1.4)*yscale(.63)*shift(-path_middle(itcont6))*itcont6;
filldraw(ugscont1, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Definition", a = (-95, -163));
picture sub16;
drawsegment(sub16, (0,-.3), (4,-.3));
label(sub16, "\Huge \(0\)", (0,-.3), align = 2.7*S);
label(sub16, "\Huge \(1\)", (4,-.3), align = 2.7*S);
draw(sub16, (6,-.3)--(8,-.3), p = p, arrow = ArcArrow(SimpleHead), L = Label("\Huge \(f\)", position = MidPoint, align = 3*S));
pair a = (10,-1);
pair b = (16,1);
draw(sub16, a--b, textbubblepen+blue);
draw(sub16, a{E+S}..(13, -.8){N}..b, p);
draw(sub16, a{E}..{N+1.2*E}b, p);
draw(sub16, a..(11.5, 1.1){E+.5*N}..{E+.1*N}b, p);
dot(sub16, Label("\Huge \(a\)", align = 2.2*S+W), a, dotpen);
dot(sub16, Label("\Huge \(b\)", align = 2.2*S+E), b, dotpen);
add(shift(-122, -171)*scale(5)*sub16);

path ugscont2 = shift(12, -163)*yscale(.81)*rotate(-5)*shift(-path_middle(iscont4))*iscont4;
filldraw(ugscont2, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Lemma", a = (0, -150));
picture sub17;
pair a1 = (0,0), b1 = (8, 1);
path sp = a1{N+.7*E}..(.35*b1)..{N+.2*E}b1;
real t = .6;
pair ft = point(sp, t*length(sp));
draw(sub17, a1--b1--ft--cycle, p);
draw(sub17, sp, textbubblepen+blue);
dot(sub17, Label("\Huge \(a\)", align = 2.5*S+W), a1, dotpen);
dot(sub17, Label("\Huge \(b\)", align = 1.5*S+2*E), b1, dotpen);
dot(sub17, Label("\Huge \(f(t)\)", align = 3*S+1*W), ft, dotpen);
add(shift(-10, -163)*scale(5.5)*sub17);
label("\Huge \(d(a, b) = d(a, f(t)) + d(f(t), b)\)", (20, -175));

path ugscont3 = shift(95, -162)*xscale(1.05)*yscale(.9)*rotate(3)*reflect((0,0), (1,0))*shift(-path_middle(iscont1))*iscont1;
filldraw(ugscont3, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Lemma", a = (95, -148));
label("\Huge \textit{In a UGS all metric segments are free polytopes.}", (95, -155));
picture sub18;
drawsegment(sub18, (0,0), (6,0));
draw(sub18, (-2,0)..(3,1.7)..(8,0)..(3,-1.6)..cycle, dashedpen);
sc = 5;
rad = 1.5/sc;
draw(sub18, circle((0,0), rad), red+2);
draw(sub18, circle((6,0), rad), red+2);
add(shift(79, -169)*scale(sc)*sub18);


// ? filling the induced structure

path istcont1 = shift(100, 115)*scale(1)*reflect((0,0), (0,1))*shift(-path_middle(itcont1))*itcont1;
filldraw(istcont1, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Topology", a = (95, 135));
ystep = -5;
xpos = 69;
ypos = 126;
label("\Huge (1) \textit{The Polytope Intersection Lemma}", (xpos, ypos), align = E);
label("\Huge (2) \textit{The Polytope Union Lemma}", (xpos, ypos+ystep), align = E);
label("\Huge (3) \textit{Finite dimension}", (xpos, ypos+2*ystep), align = E);
label("\Huge (4) \textit{Freedom}", (xpos, ypos+3*ystep), align = E);

picture sub19;
pair tq1 = (0,0);
pair tq2 = (3,1);
pair tq3 = (6,-1);
pair tq4 = (5,-4);
pair tq5 = (1,-4);
sc = 4.5;
rad = .5/sc;
fill(sub19, tq1--tq2--tq3--tq4--tq5--cycle, lightgrey);
draw(sub19, tq1--tq2--tq3--tq4--tq5--tq1--tq3--tq5--tq2--tq4--cycle, dashedpen);
filldraw(sub19, circle(tq1, rad), fillpen = white, drawpen = p);
filldraw(sub19, circle(tq2, rad), fillpen = white, drawpen = p);
filldraw(sub19, circle(tq3, rad), fillpen = white, drawpen = p);
filldraw(sub19, circle(tq4, rad), fillpen = white, drawpen = p);
filldraw(sub19, circle(tq5, rad), fillpen = white, drawpen = p);
add(shift(106, 113)*scale(sc)*sub19);

path istcont2 = shift(81.5,80)*rotate(-91)*scale(1.15)*shift(-path_middle(iscont3))*iscont3;
filldraw(istcont2, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "TPIL", a = (82, 97));
picture sub20;
fill(sub20, u1[0]--u1[1]--in1--u2[0]--in2--u1[3]--u1[4]--cycle, p = mediumblue);
fill(sub20, u2[1]--u2[2]--u2[3]--in2--u1[2]--in1--cycle, p = mediumyellow);
fill(sub20, in1--u1[2]--in2--u2[0]--cycle, p = mediumgreen);
drawpolytope(sub20, u1);
drawpolytope(sub20, u2);
dot(sub20, in1, dotpen);
draw(sub20, circle(in1, 1.5), 2+red);
dot(sub20, in2, dotpen);
draw(sub20, circle(in2, 1.5), 2+red);
add(shift(124, 0)*rotate(-90)*sub20);

path istcont3 = shift(99, 46.5)*scale(1.4)*xscale(.75)*rotate(30)*shift(-path_middle(iscont5))*iscont5;
filldraw(istcont3, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Local", a = (124.5, 74));
addtext(s = "Convexity", a = (122.5, 67));




smooth s = smooth(
    contour = reflect((0,0), (0,1))*samplesmooth(2).contour,
    holes = new hole[]{
        hole(
            contour = samplesmooth(2).holes[0].contour,
            shift = (.1,-.6),
            scale = 1.4,
            rotate = -70,
            sections = new real[][]{
                new real[] {dummynumber, dummynumber, 290, 10, .65, 200}
            }
        )
    },
    subsets = new subset[]{
        subset(
            contour = yscale(1.3)*shift(-path_middle(iscont3))*iscont3,
            scale = .08,
            rotate = 70,
            shift = (2,-.7),
            label = "\Huge \(U\)",
            labeldir = W+.2*N
        )
    }
);
s.move(shift = (97, 40), rotate = 50, scale = 5.3);
s.set_label("\Huge \(M\)", labeldir = (W+.4*S));
draw(s, contourpen = textbubblepen, sectionpen = p, viewdir = dir(45), subsetpen = mediumgreen+cyan);





filldrawpolytope(new pair[]{(100, 43), (107, 50), (112, 44)});


path istcont4 = shift(98.5, -4)*rotate(20)*yscale(.85)*scale(1.05)*shift(-path_middle(iscont4))*iscont4;
filldraw(istcont4, fillpen = white, drawpen = boxtitlebubblepen);
plainaddtext(s = "\Huge \(E^2 \ne B^2\)", a = (91, -21), offset = 2.3);
picture sub21;
draw(sub21, (-3,3)--(3,3)--(3,-3)--(-3,-3)--cycle, textbubblepen);
filldrawpolytope(sub21, scale(.95)*sequence(new pair(int i){return (cos(i/5*2*pi), sin(i/5*2*pi));}, 5));
drawpolytope(sub21, new pair[]{(-2.7, -.2), (1, 2.3), (2.1, -2.5)}, drawpen = dashedpen);
add(shift(111, 0)*scale(5)*sub21);
picture sub22;
filldraw(sub22, scale(1.8)*unitcircle, drawpen = dashedpen, fillpen = lightgrey);
filldrawpolytope(sub22, scale(1.55)*sequence(new pair(int i){return (cos(i/5*2*pi), sin(i/5*2*pi));}, 5));
add(shift(73, -8)*scale(5)*sub22);
draw((111,0)--(73,-8), p = p+red, arrow = ArcArrows(SimpleHead), margin = Margin(38, 25));

path istcont5 = shift(97, -52)*scale(1.03)*yscale(.97)*shift(-path_middle(itcont6))*itcont6;
filldraw(istcont5, fillpen = white, drawpen = boxtitlebubblepen);
addtext(s = "Riemannian", a = (96.5, -35));
picture sub23;
smooth s2 = samplesmooth(3);
s2.move(scale = 3.6, rotate = -90);
draw(sub23, s2, viewdir = dir(-45), contourpen = p, sectionpen = linewidth(.4mm));
pair clippair = (15,2);
path clippath = circle(clippair, 4);
draw(sub23, clippath, p+blue);
picture sub24;
sub24.add(sub23);
clip(sub24, clippath);
transform clipsh = shift(-6.5, -3)*scale(3);
pair clippair2 = clipsh*clippair;
pair x4 = (13.5,0);
pair y4 = (17.4,3.7);
draw(sub24, curved_path(x4, y4, curve = .1), p+blue);
dot(sub24, Label("\Huge \(x\)", align = 1.5*E+1.5*S), x4, dotpen);
dot(sub24, Label("\Huge \(y\)", align = 1.5*E+S), y4, dotpen);
sub23.add(clipsh*sub24);
draw(sub23, clipsh*clippath, p+blue);
add(shift(81.5,-56)*sub23);
addtext(s = "LUGS", a = (101, -67));

// ? filling references

path refcont = shift(97, -114.5)*xscale(.85)*yscale(.65)*shift(-path_middle(refbubble))*refbubble;
filldraw(refcont, fillpen = white, drawpen = boxtitlebubblepen);
xpos = 65;
ypos = -105.5;
real xstep = 4;
real ystep1 = -3;
real ystep2 = -5;

label("\LARGE (1) M.L.J. van de Vel, \textit{Theory of convex structures},", (xpos, ypos), align = E);
label("\LARGE North-Holland mathematical library, 1993;", (xpos+xstep, ypos+ystep1), align = E);
label("\LARGE (2) D. C. Kay, E. W. Womble, \textit{Axiomatic convexity theory},", (xpos, ypos+ystep1+ystep2), align = E);
label("\LARGE Pacific Journal of Mathematics, 1971;", (xpos+xstep, ypos+2*ystep1+ystep2), align = E);
label("\LARGE (3) D. Gromoll, W. Klingenberg and W. Meyer, \textit{Riemannsche}", (xpos, ypos+2*ystep1+2*ystep2), align = E);
label("\LARGE \textit{Geometrie im Grossen}, Springer-Verlag, Berlin, 1968.", (xpos+xstep, ypos+3*ystep1+2*ystep2), align = E);

// ? Final step!

// addtext(s = "\textit{Roman Maksimovich}", a = (c3-10, (bodyd-h2)/2), inches = .4);
addtext(s = "\textit{Made with} \texttt{Asypmtote}", a = (c4-40, (bodyd-h2)/2), inches = .4);
