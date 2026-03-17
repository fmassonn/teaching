Conversion en PDF: 

/opt/anaconda3/bin/python3 -m jupyter nbconvert --to pdf Chap07.ipynb --output Chap07

et voir les erreurs dans la console

erreurs typiques:

il faut aller à la ligne après $$: 

$$
x=3
$$

et ne pas laisser de ligne blanche dans les tableaux/accolades etc.

dans le inline text, ne pas laisser d'espace (pas "$x = Y + 4 $ mais $x = Y + 4$)
