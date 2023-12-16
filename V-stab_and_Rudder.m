ssalpha = 20/57.3;
syms sigma deflcw  % Separate declarations for each symbolic variable
Ixxb = 3;
Izzb = 4.24;
Ixzb = 0.04;
Ixxw = Ixxb*((cos(ssalpha))^2)+Izzb*((sin(ssalpha))^2)-Ixzb*sin(2*ssalpha);
Izzw = Ixxb*((sin(ssalpha))^2)+Izzb*((cos(ssalpha))^2)+Ixzb*sin(2*ssalpha);
Ixzw = Ixxb*0.5*sin(2*ssalpha)-Izzb*0.5*sin(2*ssalpha)+Ixzb*cos(2*ssalpha);
RSR = 1.4;
NSR = ((Ixxw*Izzw - Ixzw^2)/Ixxw)*RSR;
fprintf('NSR: %.6f\n', NSR);

T = 0;

b = 3;
Svmin = 100;

defr = -100;
bv = 0.29;
bv1 = 100;
Cvr = 0.34;
 % Initialize tcv before the outer loop
Cvt1 = 1000;
Crud1 = 100;
Srud1 = 100;
tcv1=100;
Cnb1=100;
Cndr1=100;
Cydr1=100;
Cyb1=100;
eqn3=100;
eqn4=100;
eqn5=100;
eqn6=100;
eqn7=100;
eqn8=100;
eqn9=100;
eqn10=100;
eqn11=100;
sigma1=100;
deflcw1=100;
defr1=100;

while true
    bv = bv + 0.01;
    tcv = 0.5;
    if bv>0.5
        break;
    end 
    while true
        tcv = tcv + 0.01;
        Cvt = tcv * Cvr;
        %fprintf('tcv: %.6f\n', tcv);
        Sv = (Cvr + Cvt) * 0.5 * bv;
        CV = Sv / bv;
        SVe = Sv - 0.3 * bv * CV;
        %fprintf('Sve: %.6f\n', SVe);
        X = 1.35;
        MACv = (2/3) * Cvr * ((1 + tcv + tcv^2) / (1 + tcv));
        Xcg = 0.427;
        lv = X + (Cvr - 0.75 * MACv) - Xcg;
        %fprintf('lv: %.6f\n', lv*1);
        %fprintf('S: %.6f\n', 1.81);
        %fprintf('b: %.6f\n', b);
        clav = 0.10391 * 57.3;
        ARv = 2 * (bv^2) / Sv;
        ARw = 4.76;
        Clav = clav / (1 + ((clav) / (3.14 * ARv)));
        Vve= double((lv * SVe) / (b * 1.81));
    
        % Extract the numerical value from the solution struct
        
    
        %fprintf('Vve: %.6f\n', Vve);


        %fprintf('Vve: %.6f\n', Vve);
        %Crud = Cvt - 0.09;
        Crud=0.136;
        Srud = Crud * bv;
        tau = 1.129 * ((Srud / Sv)^0.4044) - 0.1772;
        n = 0.95;
        Vs = 11;
        Cndr = -Clav * Vve * n * tau * 1;
        %fprintf('Cndr: %.6f\n', Cndr);
        defr = (2 * NSR) / (0.768 * (Vs^2) * 1.81 * b * Cndr);
        %fprintf('bv: %.6f\n', bv);
        %fprintf('Cvt: %.6f\n', Cvt);
        %fprintf('defr (deg): %.6f\n', defr * 57.3);
        Vv = (Sv*lv)/(1.81*b);
        Crw = 0.645;
        Ctw = 0.4;
        sweep = atan(3*(Crw-Ctw)/(4*b));
        nv = 0.724 + 3.06*((Sv/1.81)/(1+cos(sweep)));
        Cnb = nv*Vv*Clav;
        %fprintf('Cnb: %.6f\n', Cnb);
        Vw = 9;
        VT = sqrt(Vw^2 + Vs^2);
        %beta = atan(Vw/Vs);
        beta=30/57.3;
        %fprintf('beta (deg): %.6f\n', beta*57.3);
        Cdy = 0.5;
        Ss = 0.12 + Sv;
        Fw = 0.5*1.225*(Vw^2)*Ss*Cdy;
        xf = 0.419;
        xv = X + ((2*Cvt*(Cvr-Cvt)+(Cvt^2)+(Cvr-Cvt)*Cvr+Cvt*Cvr+Cvr^2)/(3*(Cvr+Cvt)));
        xca = (xf*0.12 + xv*Sv)/(xf+Sv);
        dc = Xcg - xca;
        Cl = 0.16;
        %Cyb = (Cl^2)*((6*tan(sweep)*sin(sweep))/(3.14*ARw*(ARw+4*cos(sweep))));
        Cyb = -n*(Sv/1.81)*Clav*1;
        Cydr = Clav*(Sv/1.81)*tau;
        Cndr = -Clav * Vv * n * tau * 1;
        % Update the symbolic solver part with vpasolve
        eqn1 = 0.5*1.225*(VT^2)*1.81*b*(2*Cnb*(beta-sigma)+Cndr*deflcw)+Fw*dc*cos(sigma) == 0;
        eqn2 = 0.5*1.225*(Vw^2)*Ss*Cdy-0.5*1.225*(VT^2)*1.81*(Cyb*(beta-sigma)+Cydr*deflcw) -T*sin(sigma)== 0;

        %eqn1 = Cyb*(beta - sigma) + Cydr*deflcw == 0;
        %eqn2 = Cnb*(beta - sigma) + Cndr*deflcw == 0;
        
        S = vpasolve([eqn1, eqn2], [sigma, deflcw]);
        %fprintf('deflcw: %.6f\n', S.deflcw);
        %fprintf('sigma (deg): %.6f\n', S.sigma*57.3);
        

        
        if tcv <= 0.8 && Cnb >= 0.05 && S.deflcw <= 2*30/57.3 && defr >= -2*30/57.3 && Cvt-Crud>=0.09
            %fprintf('Sv: %.6f\n', Sv);
            if Sv<Svmin 
                Svmin=Sv;
                Crud1 = Crud;
                Srud1 = Srud;
                Cvt1 = Cvt;
                bv1 = bv;
                tcv1=tcv;

                Cnb1=Cnb;
                Cydr1=Cydr;
                Cndr1=Cndr;
                Cyb1=Cyb;
                
                sigma1=S.sigma;
                deflcw1=S.deflcw;
                defr1=defr;
            end
        end

        if tcv > 0.8
            break;  % Break out of the inner loop when tcv exceeds 0.8
        end
          % Update tcv within the inner loop
    end
end
fprintf('tcv: %.6f\n', tcv1);
fprintf('def spin recovery (deg): %.6f\n', defr1 * 57.3);

fprintf('beta (deg): %.6f\n', beta*57.3);
fprintf('defl cross wind (deg): %.6f\n', deflcw1*57.3);
fprintf('sigma (deg): %.6f\n', sigma1*57.3);
fprintf('Cvt: %.6f\n', Cvt1);
fprintf('Crud: %.6f\n', Crud1);
fprintf('bv: %.6f\n', bv1);
fprintf('Sv: %.6f\n', Svmin);

fprintf('Srud/Sv: %.6f\n', Srud1/Svmin);

fprintf('Cnb: %.6f\n', Cnb1);
fprintf('Cyb: %.6f\n', Cyb1);
fprintf('Cndr: %.6f\n', Cndr1);
fprintf('Cydr: %.6f\n', Cydr1);



fprintf('T*sin(sigma): %.6f\n', eqn11);



