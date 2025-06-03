# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 20:54:20 2024

@author: vador
"""

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

# Stockage de l'image
image = Image.open('pillars_of_creation.jpg')
plt.figure(figsize=(8, 8))
plt.imshow(image)
plt.show()


def svd_image(image, order):
    # Stockage de la matrice compressée
    compressed = np.zeros(image.shape)
    
    # V : ensemble de vecteurs de base orthonormés de Kn, vecteur singulier à droite
    # U : ensemble de vecteurs de base orthonormés de Km, vecteur singulier à gauche
    # Σ : valeurs singulières i.e. racines des valeurs propres (triée par ordre décroissant)
    U, S, V = np.linalg.svd(image)
    
    # Reconstruction de la matrice compressée
    for i in range(order):
        Ui = U[:, i].reshape(-1, 1)
        Vi = V[i, :].reshape(1, -1)
        Si = S[i]
        compressed += (Ui * Si * Vi)
    
    return compressed


# Conversion de l'image en degrés de gris uniquement
gray_image = image.convert('L')
plt.figure(figsize=(8, 8))
plt.imshow(gray_image, cmap='gray')
plt.show()

# Conversion en numpy array
gray_image = np.array(gray_image)

orders = [1, 5, 10, 20, 50, 100, 200, 300, 400, 500]
erreur = []

plt.figure(figsize=(25, 10))

# Boucle sur les différents ordres
for i in range(len(orders)):
    order = orders[i]
    gray_comp = svd_image(gray_image, order)
    
    # Arrondi à l'unité près (0,255)
    gray_comp = np.around(gray_comp).astype(int)
    
    # Calcul de l'erreur
    erreur.append(np.sqrt(np.mean((gray_image - gray_comp)**2)))
    
    # Plot
    plt.subplot(2, 5, i + 1)
    plt.title("Ordre = {} | Erreur = {:.2f}".format(order, erreur[i]))
    plt.axis('off')
    plt.imshow(gray_comp, cmap='gray')

plt.suptitle('Compression aux différents ordres')
plt.show()

plt.figure(figsize=(15, 4))
plt.scatter(orders, erreur, label="Images compressées")
plt.plot(orders, erreur, c='r', ls='-', label="fit")
plt.title("Diminution de l'erreur avec l'ordre")
plt.xlabel("Ordre de Compression")
plt.ylabel("Erreur")
plt.grid('on')
plt.legend()
plt.show()



original_size = 575 * 600 * 1
orders = np.arange(1, 575 + 1, 1)
compression_size = (575 + 1 + 600) * orders

compression_ratio = compression_size / original_size

plt.figure(figsize=(15, 4))

plt.plot(orders, compression_ratio, label="Courbe de ratio")
plt.plot(orders, np.ones((575)), c='r', ls='--', label="limite")
plt.fill_between(orders, compression_ratio, 1, where=(orders < 288),
                 color='green', alpha=0.1, label="Gain")             
plt.fill_between(orders, compression_ratio, 1, where=(orders >= 288),
                 color='red', alpha=0.1, label="Perte")   
    
plt.title("Evolution du ratio de compression en fonction de l'ordre retenu")
plt.xlabel("Ordre de Compression")
plt.ylabel("Ratio de Compression")
plt.grid('on')
plt.legend()
plt.show()


# Séparation des RGB
red_image = np.array(image)[:, :, 0]
green_image = np.array(image)[:, :, 1]
blue_image = np.array(image)[:, :, 2]

# Compression
order = 500
red_comp = svd_image(red_image, order)
green_comp = svd_image(green_image, order)
blue_comp = svd_image(blue_image, order)

# Recombinaison des images colorées
color_comp = np.zeros((np.array(image).shape[0], np.array(image).shape[1], 3))
color_comp[:, :, 0] = red_comp
color_comp[:, :, 1] = green_comp
color_comp[:, :, 2] = blue_comp
color_comp = np.around(color_comp).astype(int)

# Affichage des 3 couleurs
plt.figure(figsize=(24, 8))
plt.subplot(141)
plt.imshow(red_comp, cmap='Reds_r')
plt.title("Red Channel")
plt.subplot(142)
plt.imshow(green_comp, cmap='Greens_r')
plt.title("Green Channel")
plt.subplot(143)
plt.imshow(blue_comp, cmap='Blues_r')
plt.title("Blue Channel")
plt.subplot(144)
plt.imshow(color_comp)
plt.title("Combined RGB Channels")
plt.show()



plt.figure(figsize=(25, 10))

orders = [1, 5, 10, 20, 50, 100, 200, 300, 400, 500]

for i in range(len(orders)):
    order = orders[i]
    red_comp = svd_image(red_image, order)
    green_comp = svd_image(green_image, order)
    blue_comp = svd_image(blue_image, order)
    
    # Combinaison des images
    color_comp = np.zeros((np.array(image).shape[0], np.array(image).shape[1], 3))
    color_comp[:, :, 0] = red_comp
    color_comp[:, :, 1] = green_comp
    color_comp[:, :, 2] = blue_comp
    color_comp = np.around(color_comp).astype(int)
    
    # Plot
    plt.subplot(2, 5, i + 1)
    plt.title("Ordre = {}".format(order))
    plt.axis('off')
    plt.imshow(color_comp)

plt.suptitle('Compression aux différents ordres')
plt.show()
