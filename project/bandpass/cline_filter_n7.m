% Gleichungen für die Berechnung eines Coupled Line Bandpass Filter
% Die Gleichungen gelten für Microstrip Leitungen
%
% Berechnung eines Coupled Line Kopplers
% --------------------------------------
% Die Berechnung ist als erster Wurf für eine Optimierung zu betrachten.
% Es werden die Gleichungen und das Vorgehen nach LEE (2004, 216-217)
% angewendet.
%
% Jom, 17.4.2012


% Feldkonstanten
mu_0 = 4*pi*1e-7;
e_0 = 8.854e-12;

% Variablen des Printmaterials (z.B. RO_4350 B)
T = 17e-6;      % Dicke der Kupferkaschierung (m)
sigma = 5.88e7; % Leitfähigkeit des Leitermaterials (Sm)

e_r = 3.66;      % rel. Permittivität des Dielektrikums
mu_r = 1;        % rel. Permeabilität des Dielektrikums
H = 1.524e-3;   % Dicke des Dielektrikums (m)
tanD = 0.004;   % Verlustwinkel tan_delta des Materials

% Variablen für den Entwurf
f_0 = 2.449e9;    % Entwurfsfrequenz 1.65 GHz
z_0 = 50;       % Systemimpedanz 50 Ohm (Bezugsimpedanz)
B = 80e6;

% Filterkoeffizienten
% Tschebyscheff 0.5 dB 5. Ordnung
G = [1.737 1.258 2.638 1.344 2.638 1.258 1.737];


% Berechnung der Admittanzinverter J
J =0;
J(1) = 1/z_0*sqrt(pi*B/2/f_0/G(1));
for i=2:length(G),
    J(i)=1/2/z_0*pi*B/f_0/sqrt(G(i-1)*G(i));
end
J(length(G)+1)=1/z_0*sqrt(pi*B/2/f_0/G(length(G))*1);

% Berechnung der Z_even und Z_odd Impedanzen der einzelnen Koppler
CLine = 0;
for i = 1:length(J),
    CLine(i,1) = z_0*(1+J(i)*z_0+(J(i)*z_0)^2); % Z_even_i
    f = @(x)sqrt(mu_0*mu_r/e_0/e_r)/x*(1+1.735*e_r^(-0.0724)*(x)^(-0.836))^(-1)-CLine(i,1)/2;
    we_h = fzero(f,10);
    
    CLine(i,2) = z_0*(1-J(i)*z_0+(J(i)*z_0)^2); % Z_odd_i
    f = @(x)sqrt(mu_0*mu_r/e_0/e_r)/x*(1+1.735*e_r^(-0.0724)*(x)^(-0.836))^(-1)-CLine(i,2)/2;
    wo_h = fzero(f,10);
    
    % BeRechnung der Startwerte für die Iteration 
    s_h = 2/pi*acosh((cosh(pi*wo_h/2)+cosh(pi*we_h/2)-2)/(cosh(pi*wo_h/2)-cosh(pi*we_h/2)))
    w_h = we_h/2

    error = 1e-6; %
    eps1 = 0.1;
    eps2 = 0.1;
    while (eps1 > error)|(eps2 >error),
        g = cosh(pi*s_h/2);
        d = cosh(pi*(w_h + s_h/2));
    
        we_h2 = 2/pi*acosh((2*d-g+1)/(g+1));
    
        if e_r < 6
            wo_h2 = 2/pi*acosh((2*d-g-1)/(g-1)) + 4/pi/(1+e_r/2)*acosh(1+2*w_h/s_h);
        else
            wo_h2 = 2/pi*acosh((2*d-g-1)/(g-1)) + 1/pi*acosh(1+2*w_h/s_h);
        end
        eps1 = ((we_h2 + wo_h2) - (we_h + wo_h)) / (we_h+wo_h);
        w_h = w_h - 1/2*eps1*w_h;
    
        eps2 = ((we_h2 - wo_h2) - (we_h - wo_h)) / (we_h-wo_h);
        s_h = s_h + 1/2*eps2*s_h;
    end
    % Resultatausgabe
    CLine(i,3) = s_h*H*1000;  % S Wert
    CLine(i,4) = w_h*H*1000;  % H Wert
end
% Resultatausgabe
% Ausgabeformat: Z_even , Z_odd, S, W
CLine 



