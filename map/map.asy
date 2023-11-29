settings.outformat = "pdf";
import roundedpath;
texpreamble("\input{Math}");
texpreamble("\usepackage{color}");

path connection_arc (pair a = (0,0), pair b = (1, 0), real curve = 0)
{
    pair m = ((a+b)/2)+(unit(rotate(-90) * (b-a)) * curve);

    return a .. m .. b;
}

path connection_arrow (pair a = (0,0), pair b = (1, 0), pair start = (1, 0), pair end = (1, 0))
{
    return a{start} .. {end}b;
}

size(0, 27.5cm);

pair tangencePoint_1 = (-8, 6);
pair tangencePoint_2 = (10, 4);
pair tangencePoint_3 = (7.5, -6.5);
pair tangencePoint_4 = (-6, -6.5);

pair tangenceVector_1 = (1,1);
pair tangenceVector_2 = (3,-5);
pair tangenceVector_3 = (-3,-2);
pair tangenceVector_4 = (-3,3);

// % ! homotopy line

pair homotopyStart = tangencePoint_1 + (-15,-5);
pair homotopyEnd = tangencePoint_1 + (-4,13);

real homotopyNoteHeight = 1.3;
real homotopyNoteWidth = 7;

path homotopyBox = roundedpath(homotopyEnd+(homotopyNoteWidth/2, 0) -- homotopyEnd+(-homotopyNoteWidth/2, 0) -- homotopyEnd+(-homotopyNoteWidth/2, -homotopyNoteHeight) -- homotopyEnd+(homotopyNoteWidth/2, -homotopyNoteHeight) -- cycle, homotopyNoteWidth);

label("{\Large \textbf{H\ O\ M\ O\ T\ O\ P\ Y}}", homotopyEnd+(0, -homotopyNoteHeight/2));

pair homotopyVerticalPoint = tangencePoint_1+(6, 9);

path homotopy = homotopyStart .. {tangenceVector_1}tangencePoint_1 .. {(0,1)}(homotopyVerticalPoint) .. {(-1,0)}(homotopyEnd+(homotopyNoteWidth/2, 0)) -- homotopyEnd+(-homotopyNoteWidth/2, 0){(-1,0)} .. homotopyEnd+(-homotopyNoteWidth/2 - 3, -1);

draw(homotopy);
draw(homotopyBox);

// % ! topology line

pair topologyStart = tangencePoint_2 + (-2,13);
pair topologyEnd = tangencePoint_2 + (12,-5);

real topologyNoteHeight = 1.3;
real topologyNoteWidth = 7;

path topologyBox = roundedpath(topologyEnd+(topologyNoteWidth/2, 0) -- topologyEnd+(-topologyNoteWidth/2, 0) -- topologyEnd+(-topologyNoteWidth/2, topologyNoteHeight) -- topologyEnd+(topologyNoteWidth/2, topologyNoteHeight) -- cycle, topologyNoteWidth);

label("{\Large \textbf{T\ O\ P\ O\ L\ O\ G\ Y}}", topologyEnd+(0, topologyNoteHeight/2));

pair topologyVerticalPoint = tangencePoint_2+(-3, 9);

path topology = topologyStart .. {(0,-1)}topologyVerticalPoint .. {tangenceVector_2}tangencePoint_2 .. {(1,0)}(topologyEnd+(-topologyNoteWidth/2, 0)) -- topologyEnd+(topologyNoteWidth/2, 0){(1,0)} .. topologyEnd+(topologyNoteWidth/2 + 1, .2);

draw(topology);
draw(topologyBox);

// % ! analysis line

pair analysisStart = tangencePoint_3 + (13,0);
pair analysisEnd = tangencePoint_3 + (-6.5,-13);

real analysisNoteHeight = 1.3;
real analysisNoteWidth = 6.3;

path analysisBox = roundedpath(analysisStart+(analysisNoteWidth/2, 0) -- analysisStart+(-analysisNoteWidth/2, 0) -- analysisStart+(-analysisNoteWidth/2, -analysisNoteHeight) -- analysisStart+(analysisNoteWidth/2, -analysisNoteHeight) -- cycle, analysisNoteWidth);

label("{\Large \textbf{A\ N\ A\ L\ Y\ S\ I\ S}}", analysisStart+(0, -analysisNoteHeight/2));

path analysis = analysisStart+(analysisNoteWidth/2+3, -1) .. {(-1,0)}analysisStart+(analysisNoteWidth/2, 0) -- (analysisStart+(-analysisNoteWidth/2, 0)){(-1,0)} .. {tangenceVector_3}tangencePoint_3  .. analysisEnd;

draw(analysis);
draw(analysisBox);

// % ! order line

pair orderStart = tangencePoint_4 + (3,-12);
pair orderEnd = tangencePoint_4 + (-10, 2.5);

real orderNoteHeight = 1.3;
real orderNoteWidth = 4.5;

path orderBox = roundedpath(orderEnd+(orderNoteWidth/2, 0) -- orderEnd+(-orderNoteWidth/2, 0) -- orderEnd+(-orderNoteWidth/2, -orderNoteHeight) -- orderEnd+(orderNoteWidth/2, -orderNoteHeight) -- cycle, orderNoteWidth);

label("{\Large \textbf{O\ R\ D\ E\ R}}", orderEnd+(0, -orderNoteHeight/2));

path order = orderStart .. {tangenceVector_4}tangencePoint_4 .. {(-1,0)}(orderEnd+(orderNoteWidth/2, 0)) -- orderEnd+(-orderNoteWidth/2, 0){(-1,0)} .. orderEnd+(-orderNoteWidth/2-3, -2);

draw(order);
draw(orderBox);

// % ! convex line

pair upMaxPoint = (4, 11);
pair rightMaxPoint = (13, -3.5);
pair downMaxPoint = (-1, -11);
pair leftMaxPoint = (-12, 0);

pair upMaxVector = (1,-.3);
pair rightMaxVector = (-1,-1);
pair downMaxVector;
pair leftMaxVector;

draw(tangencePoint_1{tangenceVector_1} .. {upMaxVector}upMaxPoint .. tangencePoint_2{tangenceVector_2} .. {rightMaxVector}rightMaxPoint .. tangencePoint_3{tangenceVector_3} .. downMaxPoint .. tangencePoint_4{tangenceVector_4} .. leftMaxPoint .. cycle);

// % ! up additional line

draw(upMaxPoint{-upMaxVector} .. {(0,1)}homotopyVerticalPoint .. (3,18) .. {(0, -1)}topologyVerticalPoint .. cycle);

// % ! right additional line

draw(rightMaxPoint{-rightMaxVector} ..  topologyEnd{(1,0)} .. (27, -4.5) .. {(-1,0)}analysisStart+(analysisNoteWidth/2, 0) -- analysisStart+(-analysisNoteWidth/2, 0){(-1,0)} .. cycle);

// % ! convex objects

pair textAlign = N * 0.06;
pair catAlign = N * 0.14;

pair convexPosition = (1, 3.6);
label("\(\catname{Convex}\)", convexPosition+catAlign);

pair betweennessPosition = (1, -.5);
label("\(\catname{Betweenness}\)", betweennessPosition+catAlign);

pair freePolytopesPosition = (-7, .5);
label("\(\catname{FreePolytopes}\)", freePolytopesPosition+catAlign);

pair dimPPosition = (6.5, 4);
label("\(\mathbf{dim}\ P\)", dimPPosition+textAlign);

pair dimXPosition = (7, 1);
label("\(\mathbf{dim}\ X\)", dimXPosition+textAlign);

pair Lin0Position = (10, -2);
label("\(\catname{0-Lin}\)", Lin0Position+catAlign);

pair Lin1Position = (6, -3);
label("\(\catname{1-Lin}\)", Lin1Position+catAlign);

pair Lin2Position = (3, -6);
label("\(\catname{2-Lin}\)", Lin2Position+catAlign);

pair LindddPosition = (1, -9);
dot(LindddPosition + S*.3, black);
dot(LindddPosition + (N+E)*.2 + S*.3, black);
dot(LindddPosition - (N+E)*.2 + S*.3, black);

pair crossProdPosition = (5, 8.5);
label("\(X \times Y\)", crossProdPosition+textAlign);

pair coProdPosition = (-5, 4.5);
label("\(X \sqcup Y\)", coProdPosition+textAlign);

pair factorPosition = (-7.5, 2.5);
label("\(X / Y\)", factorPosition+textAlign);

pair autConvexPosition = (-1.8, -7);
label("\(\mathbf{Aut_{Convex}}(X)\)", autConvexPosition+textAlign);

pair homIXPosition = (-2, 7);
label("$\mathbf{Hom_{Convex}}(\mathrm{I} \times X, Y)$", homIXPosition+textAlign);

pair homXYPosition = (-5.5, -2);
label("\(\mathbf{Hom_{Convex}}(X, Y)\)", homXYPosition+textAlign);

pair dimXcrossYPosition = (20, -4);
label("\(\mathbf{dim}(X \times Y) = \mathbf{max}(\mathbf{dim}\ X,\ \mathbf{dim}\ Y)\)", dimXcrossYPosition+textAlign);

pair aPosition = (0, 16);
pair bPosition = (2.5, 12.5);
pair xPosition = (5, 15);
label("\(A\)", aPosition+textAlign);
label("\(B\)", bPosition+textAlign);
label("\(X\)", xPosition+textAlign);
draw(aPosition -- xPosition, margin = Margin(3,3), arrow = ArcArrow(SimpleHead), p = blue);
label(Label("$\mathrm{in}$", position = MidPoint, align = N), aPosition -- xPosition);
draw(bPosition -- xPosition, margin = Margin(3,3), arrow = ArcArrow(SimpleHead),  p = yellow);
label(Label("\(f\)", position = MidPoint, align = S+E), bPosition -- xPosition);
draw(bPosition -- aPosition, margin = Margin(3,3), arrow = ArcArrow(SimpleHead), L = Label("\(\exists !\)", position = MidPoint, align = S+.6*W), p = dashed);

pair fCrossGPosition = (4, 1);
label("\(f \times g\)", fCrossGPosition+textAlign);

// % ! convex arrows

draw(connection_arc(convexPosition, betweennessPosition, curve = .4), margin = Margin(3.4, 3), arrow = ArcArrow(SimpleHead));

draw(connection_arc(betweennessPosition, convexPosition, curve = .8), margin = Margin(3, 3.4), arrow = ArcArrow(SimpleHead));

draw(connection_arc(convexPosition, dimPPosition, curve = .9), margin = Margin(7.5, 3.5), arrow = ArcArrow(SimpleHead));

draw(connection_arc(dimPPosition, dimXPosition, curve = -.5), margin = Margin(2.6, 3.6), arrow = ArcArrow(SimpleHead));

draw(connection_arc(convexPosition, crossProdPosition, curve = -.5), margin = Margin(2.5, 3.5), arrow = ArcArrow(SimpleHead), p = currentpen+dashed+red);

draw(connection_arc(convexPosition, coProdPosition, curve = .9), margin = Margin(3.9, 5.5), arrow = ArcArrow(SimpleHead), p = currentpen+dashed+blue);

draw(connection_arc(convexPosition, factorPosition, curve = .8), margin = Margin(6.8, 5), arrow = ArcArrow(SimpleHead), p = currentpen+dashed+red);

draw(connection_arc(freePolytopesPosition, convexPosition, curve = -.8), margin = Margin(5.8, 7.2), arrow = ArcArrow(SimpleHead));

draw(connection_arc(convexPosition, autConvexPosition, curve = 1.4), margin = Margin(4, 3.2), arrow = ArcArrow(SimpleHead), p = currentpen+blue);

draw(connection_arc(Lin0Position, Lin1Position, curve = -.5), margin = Margin(4, 5.4), arrow = ArcArrow(SimpleHead));

draw(connection_arc(Lin1Position, Lin2Position, curve = -.8), margin = Margin(3, 5.2), arrow = ArcArrow(SimpleHead));

draw(connection_arc(Lin2Position, LindddPosition, curve = -.7),margin = Margin(3.2, 4), arrow = ArcArrow(SimpleHead));

draw(connection_arc(betweennessPosition, Lin0Position, curve = -.2),margin = Margin(12, 5.2), arrow = ArcArrow(SimpleHead), p = currentpen+blue);

draw(connection_arc(betweennessPosition, Lin1Position, curve = .9),margin = Margin(3.5, 5.2), arrow = ArcArrow(SimpleHead), p = currentpen+blue);

draw(connection_arc(betweennessPosition, Lin2Position, curve = .6),margin = Margin(3, 4), arrow = ArcArrow(SimpleHead), p = currentpen+blue);

draw(connection_arc(convexPosition, homIXPosition, curve = .4),margin = Margin(2.8, 5), arrow = ArcArrow(SimpleHead), p = currentpen+red);

draw(connection_arrow(homXYPosition, convexPosition, start = (2,3), end = (3,1)), margin = Margin(3.5, 6), arrow = ArcArrow(SimpleHead), p = currentpen+dashed+red);

draw(connection_arc(dimXPosition, dimXcrossYPosition-(4,0), curve = -1.1), margin = Margin(6, 6), arrow = ArcArrow(SimpleHead), p = currentpen+red);

draw(connection_arc(convexPosition, (1, 13.6), curve = .5), margin = Margin(2.9, 2), arrow = ArcArrow(SimpleHead), p = currentpen+blue);

draw(connection_arc(crossProdPosition, homIXPosition, curve = 1.3), margin = Margin(5, 5), arrow = ArcArrow(SimpleHead), p = currentpen+red);

draw(connection_arc(crossProdPosition, fCrossGPosition, curve = .5), arrow = ArcArrow(SimpleHead), margin = Margin(3.3, 4), p = currentpen+blue);

// % ! homotopy objects

pair fHomGPosition = (-7, 14);
label("\(f \sim g\)", fHomGPosition+textAlign);

pair hConvexPosition = (-9.5, 9.5);
label("$\catname{hConvex}$", hConvexPosition+catAlign);

pair hTopPosition = (-13, 15);
label("\(\catname{hTop}\)", hTopPosition+catAlign);

pair pi1XPosition = (-14, 5);
label("\(\pi_1(X)\)", pi1XPosition+textAlign);

pair piNXPosition = (-15, 10.5);
label("\(\pi_n(X)\)", piNXPosition+textAlign);

pair sesPosition = (-18, 13);
label("\(\catname{SES}\)", sesPosition+catAlign);

pair lesPosition = (-19, 6);
label("\(\catname{LES}\)", lesPosition+catAlign);

// % ! homotopy arrows

draw(connection_arc(homIXPosition, fHomGPosition, curve = 1), margin = Margin(3, 4.5), arrow = ArcArrow(SimpleHead), p = currentpen+green);

draw(connection_arc(fHomGPosition, hConvexPosition, curve = -.7), margin = Margin(3, 5), arrow = ArcArrow(SimpleHead));

draw(connection_arc(hConvexPosition, hTopPosition, curve = .8), margin = Margin(3, 5.3), arrow = ArcArrow(SimpleHead), p = currentpen+dashed+blue);

draw(connection_arrow(homXYPosition+(-2,0), hConvexPosition, start = (-1,0), end = (1,2)), margin = Margin(3, 3), arrow = ArcArrow(SimpleHead), p = currentpen+blue);

draw(connection_arc(coProdPosition, pi1XPosition, curve = .8), margin = Margin(5, 5.5), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

draw(connection_arc(factorPosition, pi1XPosition, curve = -.9), margin = Margin(5, 5), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

draw(connection_arc(pi1XPosition, piNXPosition, curve = -1), margin = Margin(3.5, 3.2), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

draw(connection_arc(piNXPosition, lesPosition, curve = .9), margin = Margin(5.2, 4), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

draw(connection_arc(piNXPosition, sesPosition, curve = -1.1), margin = Margin(5.5, 4), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

draw(connection_arrow(pi1XPosition, sesPosition, start = (3,3), end = (-1,-.5)), margin = Margin(4, 5.3), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

// % ! topology objects

pair topPosition = (12, 10);
label("\(\catname{Top}\)", topPosition+catAlign);

pair t0Position = (18, 3);
label("\(\mathbf{T_0}\)", t0Position+textAlign);

pair t1Position = (13.5, 4.5);
label("\(\mathbf{T_1}\)", t1Position+textAlign);

pair manifoldPosition = (24, 3.5);
label("\(\catname{Manifold}\)", manifoldPosition+catAlign);

pair compactPosition = (19, 7);
label("\(\mathbf{Compact}\)", compactPosition+textAlign);

pair regularPosition = (25, 8.5);
label("\(\mathbf{Regular}\)", regularPosition+textAlign);

draw(topologyStart+(4, .5) .. (19, 10) .. (27, 12));

pair metPosition = (16,16);
label("\(\catname{Met}\)", metPosition+catAlign);

pair lengthMetPosition = (22, 13);
label("\(\catname{LengthMet}\)", lengthMetPosition+catAlign);

pair hausdorfPosition = (10,13);
label("\(\mathbf{Hausdorf}\)", hausdorfPosition+textAlign);

// % ! topology arrows

draw(connection_arc(convexPosition, topPosition+(-.3, 0), curve = -.8), margin = Margin(4, 4.2), arrow = ArcArrow(SimpleHead), p = currentpen+green);

draw(connection_arc(convexPosition, topPosition, curve = .4), margin = Margin(5.2, 5.5), arrow = ArcArrow(SimpleHead), p = currentpen+green);

draw(connection_arc(convexPosition, topPosition+(.4, 0), curve = 1.1), margin = Margin(6.6, 5.4), arrow = ArcArrow(SimpleHead), p = currentpen+green);

draw(connection_arc(t0Position, Lin0Position, curve = -.7), margin = Margin(2, 5), arrow = ArcArrows(SimpleHead), p = currentpen+blue);

draw(connection_arc(t1Position, Lin1Position, curve = .4), margin = Margin(2, 4), arrow = ArcArrows(SimpleHead), p = currentpen+blue);

draw(connection_arc(t0Position, t1Position, curve = 1), margin = Margin(3,3), arrow = ArcArrow(SimpleHead));

draw(connection_arc(manifoldPosition, topPosition, curve = -2), margin = Margin(9, 4), arrow = ArcArrow(SimpleHead));

draw(connection_arc(betweennessPosition, compactPosition, curve = 2.3), margin = Margin(13, 3), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

draw(point(connection_arc(betweennessPosition, compactPosition, curve = 2.3), 1.6775){(4,4)} .. {(1,2)}regularPosition, margin = Margin(0, 3), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

draw(connection_arc(topPosition, regularPosition, curve = .5), margin = Margin(4.5, 8), arrow = ArcArrow(SimpleHead));

draw(connection_arc(topPosition, compactPosition, curve = .5), margin = Margin(4.5, 8), arrow = ArcArrow(SimpleHead));

draw(connection_arc(metPosition, topPosition, curve = -.8), margin = Margin(3, 4), arrow = ArcArrow(SimpleHead));

draw(connection_arc(lengthMetPosition, metPosition, curve = .9), margin = Margin(4, 4.5), arrow = ArcArrow(SimpleHead));

draw(connection_arc(topPosition, hausdorfPosition+(.7,0), curve = .8), margin = Margin(3.3, 4.75), arrow = ArcArrow(SimpleHead));

// % ! analysis objects

pair smoothPosition = (23.5, -13);
label("\(\catname{Smooth}\)", smoothPosition+catAlign);

pair measurePosition = (15.5, -12);
label("\(\catname{Measure}\)", measurePosition+catAlign);

pair measure0Position = (11.5, -9);
label("\(\catname{0-Measure}\)", measure0Position+catAlign);

pair cKPosition = (8, -13);
label("\(C^{k}(X, Y)\)", cKPosition+textAlign);

pair lPPosition = (12, -15);
label("\(L^p(X)\)", lPPosition+textAlign);

pair lInfPosition = (21.5, -16);
label("\(L^{\infty}(X)\)", lInfPosition+textAlign);

pair cKRRPosition = (5, -17);
label("\(C^k(\R, \R)\)", cKRRPosition+textAlign);

pair uniformPosition = (16, -17);
label("\(\catname{Uniform}\)", uniformPosition+catAlign);

pair smoothPlusConvexPosition = (20, -10);
label("\(\catname{Smooth} + \catname{Convex}\)", smoothPlusConvexPosition+catAlign);

// % ! analysis arrows

draw(connection_arc(manifoldPosition, smoothPosition, curve = -6), margin = Margin(6,7), arrow = ArcArrow(SimpleHead));

draw(connection_arc(measure0Position, measurePosition, curve = -.8), margin = Margin(9, 3.1), arrow = ArcArrow(SimpleHead));

draw(connection_arc(measurePosition, lPPosition, curve = .9), margin = Margin(8,4), arrow = ArcArrow(SimpleHead));

draw(connection_arc(homXYPosition, cKRRPosition, curve = 3), margin = Margin(4,7.5), arrow = ArcArrow(SimpleHead), p = currentpen+blue);

draw(connection_arc(betweennessPosition+(-.1, 0), measure0Position, curve = 1), margin = Margin(3.8, 9.7), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

draw(connection_arc(cKRRPosition, cKPosition, curve = .6), margin = Margin(7, 4), arrow = ArcArrow(SimpleHead));

draw(connection_arc(lPPosition, lInfPosition, curve = -1), margin = Margin(6, 6), arrow = ArcArrow(SimpleHead));

draw(connection_arc(smoothPosition+(0, .2), smoothPlusConvexPosition, curve = -1.5), margin = Margin(8, 3.5), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

// %  ! order objects

pair preOrdPosition = (-10.5, -7.5);
label("\(\catname{PreOrd}\)", preOrdPosition+catAlign);

pair linOrdPosition = (-6, -13.5);
label("\(\catname{LinOrd}\)", linOrdPosition+catAlign);

pair netPosition = (-17, -8);
label("\(\mathbf{Net}\)", netPosition+textAlign);

// % ! order arrows

draw(connection_arc(preOrdPosition, betweennessPosition, curve = 1.5), margin = Margin(7.5, 4), arrow = ArcArrow(SimpleHead), p = currentpen+green);

draw(connection_arc(linOrdPosition, preOrdPosition, curve = 1), margin = Margin(4, 7), arrow = ArcArrow(SimpleHead));

draw(connection_arrow(Lin0Position, linOrdPosition, start = (0,-2), end = (-4,1)), margin = Margin(3, 7.5), arrow = ArcArrow(SimpleHead), p = currentpen+blue+dashed);

draw(connection_arc(netPosition, preOrdPosition, curve = 1), margin = Margin(5, 5), arrow = ArcArrow(SimpleHead));

// % ! legend

pair legendPosition = (-20, -12.5);
real legendLineBreak = 1;

draw(legendPosition+(-.5, legendLineBreak) -- legendPosition+(10, legendLineBreak) -- legendPosition+(10, -6) -- legendPosition+(-.5, -6) -- cycle);

label("solid line \ \(\longleftrightarrow\) \ defined concept", legendPosition, align = E);

label("dashed line \ \(\longleftrightarrow\) \ undefined concept", legendPosition+(0, -legendLineBreak), align = E);

label("\(\color{green}{\bullet}\) color \ \(\longleftrightarrow\) \ important/key", legendPosition + (0, -2*legendLineBreak), align = E);

label("\(\color{blue}{\bullet}\) color \ \(\longleftrightarrow\) \ interesting/relevant", legendPosition + (0, -3*legendLineBreak), align = E);

label("\(\color{red}{\bullet}\) color \ \(\longleftrightarrow\) \ difficult/troubling", legendPosition + (0, -4*legendLineBreak), align = E);

label("\(\color{black}{\bullet}\) color \ \(\longleftrightarrow\) \ known/given", legendPosition + (0, -5*legendLineBreak), align = E);
