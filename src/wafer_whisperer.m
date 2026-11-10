%% Wafer Whisperer
% Detection de defauts (particules, rayures) sur un wafer synthetique
% via filtrage FFT + notch des harmoniques du motif periodique.
% Reproductibilite assuree par un seed fixe.
clc; clear; close all;

% --- Parametres globaux ---
N = 512;        % taille image
P = 32;         % periode de la grille (taille d'un die en pixels)
seed = 42;
rng(seed);      % pour avoir les memes defauts a chaque run

%% Grille periodique (lignes de scribe)
% On construit le fond du wafer : un damier regulier qui simule la decoupe
% entre les dies. C'est cette periodicite qui va nous donner des pics nets
% dans la FFT, et qu'on va ensuite enlever.

background    = 0;      % intensite du fond
scribe_values = 1;      % intensite des lignes de scribe

img    = zeros(N,N);
idx    = 1:P:N;                     % positions des lignes scribe
dc_idx = floor(N/2) + 1;            % indice de la composante DC apres fftshift

%% Motif interne : reseau de vias dans chaque die
% Chaque die contient une grille reguliere de petits carres (vias).
% Ca ajoute des harmoniques supplementaires dans la FFT.

inner_period = 8;       % pixels entre deux vias
inner_size   = 2;       % taille d'un via en pixels
inner_value  = 0.5;     % intensite intermediaire

inner_idx = 1:inner_period:P;       % positions locales des vias dans un die
die = zeros(P,P);
for k = 1:length(inner_idx)
    debut_l = inner_idx(k);
    fin_l   = debut_l + inner_size - 1;
    for l = 1:length(inner_idx)
        debut_col = inner_idx(l);
        fin_col   = debut_col + inner_size - 1;
        die(debut_l:fin_l, debut_col:fin_col) = inner_value;
    end
end

% On tile le die sur toute l'image, puis on redessine les lignes scribe
% par dessus (elles ecrasent les vias en bord de die, c'est voulu).
img = repmat(die, N/P, N/P);
img(idx,:) = scribe_values;
img(:,idx) = scribe_values;

% FFT du wafer propre : sert de reference visuelle
F = fftshift(fft2(img));

figure;
imagesc(img); axis image; colorbar;
title('Wafer propre');

figure;
F_display = F;
F_display(dc_idx, dc_idx) = 0;      % on masque la DC pour la visu
imshow(log(1+abs(F_display)), []);
title('Spectre du wafer propre');

%% Particules
% Defauts ronds : simulent des contaminants (poussieres, residus de process).
% On garde la position dans ground_truth pour pouvoir evaluer le detecteur apres.

particle_radius_min = 5;
particle_radius_max = 15;
N_particles_min     = 10;
N_particles_max     = 40;
particle_value      = 0.6;

N_particles  = randi([N_particles_min, N_particles_max]);
ground_truth = false(N,N);

assert(size(img,1) == size(img,2), "The matrix should be square");
[X,Y] = meshgrid(1:N, 1:N);         % grilles de coordonnees reutilisees pour les rayures aussi

for k = 1:N_particles
    % Rayon aleatoire, centre tire dans une zone qui garde la particule
    % entierement dans l'image (evite les troncatures en bord).
    particle_radius   = randi([particle_radius_min, particle_radius_max]);
    particle_center_x = randi([particle_radius, N-particle_radius]);
    particle_center_y = randi([particle_radius, N-particle_radius]);

    % Masque disque : tous les pixels a distance <= rayon du centre.
    particle_mask = (X-particle_center_x).^2 + (Y-particle_center_y).^2 <= particle_radius.^2;

    ground_truth(particle_mask) = 1;
    img(particle_mask)          = particle_value;
end

figure;
imagesc(ground_truth); axis image; colorbar;
title("Particle ground truth");

figure;
imagesc(img); axis image; colorbar;
title("Image with particles");

%% Rayures
% Defauts allonges : segments d'orientation et de longueur aleatoires.
% On calcule le masque comme "tube" autour du segment [A,B] en projetant
% chaque pixel de l'image sur ce segment.

scratch_length_min = 50;
scratch_length_max = 200;
scratch_width_min  = 1;
scratch_width_max  = 3;
scratch_value      = 0.7;
N_scratches_min    = 7;
N_scratches_max    = 18;

N_scratches      = randi([N_scratches_min, N_scratches_max]);
scratch_mask_all = false(N,N);      % cumule toutes les rayures pour la visu

for k = 1:N_scratches
    scratch_length = randi([scratch_length_min, scratch_length_max]);
    scratch_width  = randi([scratch_width_min,  scratch_width_max]);
    scratch_angle  = 2*pi*rand;

    % Marge calculee sur la longueur max pour qu'aucune rayure ne depasse de l'image.
    scratch_margin   = floor(scratch_length_max/2);
    scratch_center_x = randi([scratch_margin, N - scratch_margin]);
    scratch_center_y = randi([scratch_margin, N - scratch_margin]);

    % Extremites A et B du segment, centrees autour de (cx,cy).
    x_a = scratch_center_x - scratch_length/2 * cos(scratch_angle);
    y_a = scratch_center_y - scratch_length/2 * sin(scratch_angle);
    x_b = scratch_center_x + scratch_length/2 * cos(scratch_angle);
    y_b = scratch_center_y + scratch_length/2 * sin(scratch_angle);

    % Projection orthogonale de chaque pixel sur le segment AB.
    % t = parametre de la projection (0 en A, 1 en B), clampe pour rester sur le segment.
    dx = x_b - x_a;
    dy = y_b - y_a;
    t  = ((X-x_a).*dx + (Y-y_a).*dy) ./ (dx^2 + dy^2);
    t  = max(0, min(1, t));

    Xproj = x_a + t.*dx;
    Yproj = y_a + t.*dy;

    % Distance euclidienne du pixel a sa projection -> masque = tube de largeur w/2.
    distance_to_scratch = hypot(X-Xproj, Y-Yproj);
    scratch_mask        = distance_to_scratch <= scratch_width/2;

    scratch_mask_all = scratch_mask_all | scratch_mask;
    ground_truth(scratch_mask) = 1;
    img(scratch_mask)          = scratch_value;
end

figure;
imagesc(img); axis image; colorbar;
title('Image avec rayures');

figure;
imagesc(scratch_mask_all); axis image; colorbar;
title('Masque cumulatif des rayures');

figure;
imagesc(ground_truth); axis image; colorbar;
title('Ground Truth - defauts artificiels');

%% Bruit
% On ajoute un bruit gaussien pour se rapprocher d'une vraie acquisition,
% puis on clippe dans [0,1] pour rester dans une plage d'intensite physique.

sigma       = 0.05;
img_noisy   = img + randn(N,N)*sigma;
img_clipped = min(1, max(0, img_noisy));

figure;
imagesc(img_clipped); axis image; colorbar;
title('Image bruitee');

% Sauvegarde du dataset pour reutilisation (training, benchmark, etc.)
save('WaferWhisperer_dataset.mat', ...
     'img','img_clipped','img_noisy','N','P','seed','sigma','ground_truth');

%% FFT de l'image degradee
% Les defauts aperiodiques brisent la periodicite -> ils vont s'etaler
% partout dans le spectre, tandis que la grille reste concentree sur ses pics.
F_defects = fftshift(fft2(img_clipped));

figure;
subplot(1,2,1);
imshow(log(1+abs(F_display)), []);
title('Clean spectrum');
subplot(1,2,2);
F_defects_display = F_defects;
F_defects_display(dc_idx, dc_idx) = 0;
imshow(log(1+abs(F_defects_display)), []);
title('Spectrum with defects + noise');

%% Filtre notch 2D
% Idee : detecter les pics intenses du spectre (= harmoniques de la grille),
% les masquer, et reconstruire par IFFT. Ce qui reste = defauts + bruit.

k_sigma      = 4;                   % seuil : moyenne + k*sigma sur le spectre
notch_radius = 3;                   % rayon de dilatation autour de chaque pic

mag_spectrum                 = abs(F_defects);
mag_spectrum(dc_idx, dc_idx) = 0;   % on ignore la DC pour le calcul du seuil

mu_spectrum    = mean(mag_spectrum, "all");
sigma_spectrum = std(mag_spectrum, 0, "all");
peaks_spectrum = mag_spectrum > (mu_spectrum + k_sigma*sigma_spectrum);

% Dilatation : on elargit chaque pic pour aussi retirer ses voisins immediats
% (fuite spectrale). Puis inversion -> notch_mask vaut 1 la ou on garde.
se         = strel('disk', notch_radius);
notch_mask = ~imdilate(peaks_spectrum, se);
notch_mask(dc_idx, dc_idx) = false; % on vire aussi la DC -> residu a moyenne nulle

F_filtered   = F_defects .* notch_mask;
img_residual = abs(ifft2(ifftshift(F_filtered)));   % module pour avoir une image >= 0

figure();
subplot(1,2,1);
imshow(notch_mask);
title('2D Notch Filter Mask');

subplot(1,2,2);
imshow(img_residual, []);
title('Residual Image (Defects + Noise)');

%% Post-traitement spatial & segmentation
% Le residu est bruite -> on lisse, on seuille statistiquement, puis on
% nettoie morphologiquement pour obtenir un masque binaire des defauts.

sigma_smooth = 1.2;     % std du filtre gaussien (px)
k_res        = 3.0;     % seuil : moyenne + k*sigma sur le residu
min_px_area  = 5;       % composantes connexes plus petites = bruit, on les enleve

% 1. Lissage gaussien : attenue le bruit haute frequence, garde les defauts
img_smooth = imgaussfilt(img_residual, sigma_smooth);

% 2. Seuillage statistique : tout pixel sortant "clairement" de la distribution
%    du fond est candidat defaut.
mu_res    = mean(img_smooth, "all");
sigma_res = std(img_smooth, 0, "all");
bw_raw    = img_smooth > (mu_res + k_res * sigma_res);

% 3. Nettoyage morphologique :
%    - imfill rebouche les trous (une particule apparait parfois en anneau)
%    - bwareaopen vire les petites composantes isolees (= bruit residuel)
bw_filled   = imfill(bw_raw, 'holes');
defect_mask = bwareaopen(bw_filled, min_px_area);

% 4. Comparaison visuelle avec la ground truth
figure('Name', 'Defect Detection vs Ground Truth');
subplot(1,3,1);
imshow(img_smooth, []);
title('Smoothed Residual');

subplot(1,3,2);
imshow(defect_mask);
title('Detected Defects (Binary)');

subplot(1,3,3);
imshow(ground_truth);
title('Ground Truth Mask');

%% Evaluation quantitative : FFT pipeline vs baseline Otsu
% On compare notre pipeline a un seuillage Otsu direct sur l'image bruitee.
% Otsu ne connait pas la structure periodique -> il va confondre la grille
% avec des defauts, c'est la ou notre methode doit gagner.

% --- A. Metriques du pipeline FFT + notch ---
TP = sum(defect_mask & ground_truth, "all");
FP = sum(defect_mask & ~ground_truth, "all");
FN = sum(~defect_mask & ground_truth, "all");

precision_fft = TP / max(TP + FP, eps);
recall_fft    = TP / max(TP + FN, eps);
f1_fft        = 2 * (precision_fft * recall_fft) / max(precision_fft + recall_fft, eps);

% --- B. Baseline : seuillage Otsu direct (sans FFT) ---
otsu_level    = graythresh(img_clipped);
baseline_mask = imbinarize(img_clipped, otsu_level);
TP_base = sum(baseline_mask & ground_truth, "all");
FP_base = sum(baseline_mask & ~ground_truth, "all");
FN_base = sum(~baseline_mask & ground_truth, "all");

precision_base = TP_base / max(TP_base + FP_base, eps);
recall_base    = TP_base / max(TP_base + FN_base, eps);
f1_base        = 2 * (precision_base * recall_base) / max(precision_base + recall_base, eps);

% Affichage du tableau comparatif
fprintf('\n PERFORMANCE COMPARISON \n');
fprintf('Method         | Precision | Recall | F1-Score\n');
fprintf('Otsu Baseline  |   %6.2f%% | %5.2f%% |  %6.2f%%\n', precision_base*100, recall_base*100, f1_base*100);
fprintf('FFT Whisperer  |   %6.2f%% | %5.2f%% |  %6.2f%%\n', precision_fft*100, recall_fft*100, f1_fft*100);

%% Overlay final
% Visualisation "rapport d'inspection" : defauts en rouge translucide
% sur l'image originale, cote a cote avec la baseline Otsu.
img_overlay = labeloverlay(img_clipped, defect_mask, 'Colormap', [1 0 0], 'Transparency', 0.3);

figure('Name', 'Final Wafer Inspection Report');
subplot(1,2,1);
imshow(baseline_mask);
title(sprintf('Otsu Baseline (F1 = %.1f%%)', f1_base*100));

subplot(1,2,2);
imshow(img_overlay);
title(sprintf('Wafer Whisperer Overlay (F1 = %.1f%%)', f1_fft*100));
