
;; return a list of map file names (which subsequently get used in jh-load-maps).
;; 
(define map-file-list
  (map (lambda (n)
	 (string-append "seed_" (number->string n) "_2FoFc.map"))
       (range 1 5))) ;; i.e. 22 to 26 inclusive.

;; return a list af molecule numbers
;; 
(define (jh-load-maps)
  (map (lambda (f) 
	 (let ((imol (handle-read-ccp4-map f 0)))
	   (set-map-colour imol 0.4 0.4 0.7)
	   imol)) 
       map-file-list))
    

(define show-next-boiling-map
  (let ((current-index 0))
    (lambda (imol-maps)

      (let ((n-maps (length imol-maps)))
	(if (> n-maps 0)
	    (begin
	      (if (= current-index (length imol-maps))
		  (set! current-index 0))
	      (display-maps (list (list-ref imol-maps current-index)))
	      (set! current-index (+ current-index 1))))))))


  
  ;; If you want to see the variance map use this (you'll need to
  ;; comment out above timeout addition if you want to see it at the
  ;; moment).
  ;; 
;; (make-variance-map imol-maps)


  (let ((time-out-handle #f)
	(imol-maps '()))
    
    (dialog-box-of-buttons 
     "Hubble-Bubble"  (cons 250 220)
     (list

      (list "Load Maps for Boiling" 
	    (lambda ()
	      (set! imol-maps (jh-load-maps))))

      (list "Boil maps" 
	    (lambda ()
	      (if (not time-out-handle)
		  (let ((handle (gtk-timeout-add 50 ;; in ms
						 (lambda ()
						   (show-next-boiling-map imol-maps)
						   #t))))
		    (set! time-out-handle handle)))))
	  
      (list "Stop Boiling" 
	    (lambda ()
	      (if time-out-handle
		  (begin
		    (gtk-timeout-remove time-out-handle)
		    (set! time-out-handle #f)))))

      (list "Contour Level Up"
	    (lambda () 
	      (for-each (lambda (imol-map)
			  (set-scroll-wheel-map imol-map)
			  (change-contour-level 1))
			imol-maps)))

      (list "Contour Level Down"
	    (lambda () 
	      (for-each (lambda (imol-map)
			  (set-scroll-wheel-map imol-map)
			  (change-contour-level 0))
			imol-maps))))

     "  Close  ")))
		  
