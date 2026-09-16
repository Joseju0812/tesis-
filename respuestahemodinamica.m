function HRF = respuestahemodinamica(t, samplerate)

% tiempo discreto
dt = 1/samplerate;
time = 0:dt:t-dt;

% parámetros doble gamma (SPM)
a1 = 6;    b1 = 1;
a2 = 16;   b2 = 1;
c  = 1/6;

% funciones gamma
h1 = (time.^(a1-1) .* exp(-time/b1)) / (gamma(a1)*b1^a1);
h2 = (time.^(a2-1) .* exp(-time/b2)) / (gamma(a2)*b2^a2);

% HRF doble gamma
HRF = h1 - c*h2;

% normalización (opcional pero recomendada)
HRF = HRF / max(HRF);

end
