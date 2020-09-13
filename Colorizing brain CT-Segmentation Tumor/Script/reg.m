fixed = dicomread('2.dcm');
moving = dicomread('4.dcm');

figure, imshowpair(moving, fixed )
title('Unregistered');

[optimizer,metric] = imregconfig('multimodal');


movingRegisteredDefault = imregister(moving, fixed, 'affine', optimizer, metric);

figure, imshowpair(movingRegisteredDefault, fixed)
title('A: Default registration')

optimizer.InitialRadius = optimizer.InitialRadius/3.5;

movingRegisteredAdjustedInitialRadius = imregister(moving, fixed, 'affine', optimizer, metric);
figure, imshowpair(movingRegisteredAdjustedInitialRadius, fixed)
title('Adjusted InitialRadius')


optimizer.MaximumIterations = 300;
movingRegisteredAdjustedInitialRadius300 = imregister(moving, fixed, 'affine', optimizer, metric);

figure, imshowpair(movingRegisteredAdjustedInitialRadius300, fixed)
title('B: Adjusted InitialRadius, MaximumIterations = 300, Adjusted InitialRadius.')