
function im=transform(refim, res)

X= res.X;

Xx = X(:,:,1);
Xy = X(:,:,2);

im= interp2(refim,  Xx, Xy );


