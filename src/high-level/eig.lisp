(in-package :magicl)


(defun hessenberg (matrix)
  "Reduce MATRIX to Hessenberg form via a sequence of unitary similarity transformations."
  (let ((m (nrows matrix))
	(n (ncols matrix))
	(vs nil)
	(H (deep-copy-tensor matrix)))
    (flet ((column-reflect! (A v k j)
	     "Reflect A[k:m,j] across V."
	     (let ((v-dot-A
		     (loop :for i :from k :below m
			   :for iv :from 0
			   :sum (* (conjugate (tref v iv 0)) (tref A i j)))))
	       (loop :for i :from k :below m
		     :for iv :from 0
		     :do (decf (tref A i j)
			       (* 2 (tref v iv 0) v-dot-A)))))
	   (row-reflect! (A v k i)
	     "Reflect A[i,k:n] across (dagger V)."
	     (let ((v*-dot-A
		     (loop :for j :from k :below n
			   :for iv :from 0
			   :sum (* (tref A i j) (tref v iv 0)))))
	       (loop :for j :from k :below n
		     :for iv :from 0
		     :do (decf (tref A i j)
			       (* 2 (conjugate (tref v iv 0)) v*-dot-A))))))
      (loop :for k :below (1- m)
	    ;; note: for Hessenberg we take V to have one less element
	    ;; than with "normal" QR.
	    :for v := (slice H (list (1+ k) k) (list m (1+ k)))
	    :for norm-v := (matrix-norm v)
	    :do (print H)
	    :unless (zerop norm-v)
	      :do (incf (tref v 0 0)
			(* norm-v (signum (tref v 0 0))))
		  (scale! v (/ (matrix-norm v)))
		  (loop :for j :from k :below n
			:do (column-reflect! H v (1+ k) j))
		  ;; for Hessenberg we are also applying similarity transformations,
		  ;; hence need to reflect rows by (dagger V).
		  (loop :for i :below m
			:do (row-reflect! H v (1+ k) i))
	    :do (push v vs))
      ;; compute Q?
      H)))
