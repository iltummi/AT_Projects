# app.py

import picamera
formats = [ 'bmp' ]
camera = picamera.PiCamera(resolution='600x490')
#camera.start_preview()
for f in formats:
    #print('Capturing format: ' + f)
    camera.capture('foto.' + f, format = f)
    #camera.stop_preview()
