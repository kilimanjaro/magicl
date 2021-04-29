(in-package :magicl)

(defun givens-entries (a b)
  "Compute the entries of the Givens matrix which rotates a vector to the x-axis.

Returns C and S such that [C -S; S C] @ [a; b] = [r; 0]."
  (cond ((zerop b)
         (values 1 0))
        ((> (abs b) (abs a))
         (let* ((tau (- (/ a b)))
                (s (/ (sqrt (+ 1 (* tau tau)))))
                (c (* s tau)))
           (values c s)))
        (t
         (let* ((tau (- (/ b a)))
                (c (/ (sqrt (+ 1 (* tau tau)))))
                (s (* c tau)))
           (values c s)))))

(defstruct givens-rotation
  "Representation of a Givens rotation between axis I and J, with C = cos(theta), S = sin(theta)."
  i j c s)


(defun left-apply-givens! (A givens &key (start-idx 0) trans)
  "Apply the Givens rotation to the matrix A, on the left (i.e. as row rotations).

Only updates columns with index >= START-IDX."
  (let ((n (ncols A))
        (i (givens-rotation-i givens))
        (j (givens-rotation-j givens))
        (c (givens-rotation-c givens))
        (s (givens-rotation-s givens)))
    (when trans
      (setf s (- s)))
    (loop :for p :from start-idx :below n
          :for t1 := (tref A i p)
          :for t2 := (tref A j p)
          :do (setf (tref A i p) (- (* c t1) (* s t2))
                    (tref A j p) (+ (* s t1) (* c t2))))
    A))

(defun right-apply-givens! (A givens &key (start-idx 0) trans)
  "Apply the Givens rotation to the matrix A, on the right (i.e. as column rotations).

Only updates rows with index >= START-IDX."
  (let ((m (nrows A))
        (i (givens-rotation-i givens))
        (j (givens-rotation-j givens))
        (c (givens-rotation-c givens))
        (s (givens-rotation-s givens)))
    (when trans
      (setf s (- s)))
    (loop :for p :from start-idx :below m
          :for t1 := (tref A p i)
          :for t2 := (tref A p j)
          :do (setf (tref A p i) (- (* c t1) (* s t2))
                    (tref A p j) (+ (* s t1) (* c t2))))
    A))

;;; TODO: - rename HESSENBERG to TRIDIAGONAL-FORM
;;;       - get Q matrix from above
;;;       - apply givens rotation to Q

(defun scale-column! (matrix idx scalar)
  (loop :for i :below (nrows matrix)
        :do (setf (tref matrix i idx)
                  (* scalar (tref matrix i idx)))
        :finally (return matrix)))

(defun hessenberg (matrix)
  "Reduce MATRIX to Hessenberg form via a sequence of unitary similarity transformations."
  (let* ((m (nrows matrix))
	 (n (ncols matrix))
	 (vs nil)
	 (H (deep-copy-tensor matrix))
         (Q (eye (list m n) :type (element-type matrix))))
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
      ;; update Q
      (loop :for k :from (- m 2) :downto 0
              :for v :in vs
              :do (loop :for j :below n
                        :do (column-reflect! Q v (1+ k) j)))
      ;; ensure H is real by absorbing phases into Q
      ;; basic idea: M = Q H Q^h = Q U^h U H U^h U Q^h where U is a diagonal
      ;; unitary matrix such that U H U^h is real tridiagonal
      (loop :with z := 1
            :for k :from 1 :below m
            :for w := (tref H k (1- k))
            :unless (zerop (imagpart w))
              :do (let* ((alpha (/ (conjugate w) (abs w)))
                         (r (realpart (* w alpha))))
                    (setf z (* z alpha))
                    ;; we "implicitly" scale row k of H by (conjugate z)
                    ;; and explicitly column k of H by z
                    (scale-column! Q k (conjugate z))
                    ;; end result of this scaling: H will be real
                    (setf (tref H k (1- k)) r
                          (tref H (1- k) k) r)))
      (values H Q))))


;; (defun hessenberg-qr! (matrix)
;;   "Compute the QR factorization of a square Hessenberg matrix, via Givens rotations.

;; Updates MATRIX, so that it is upper triangular. Returns the list of Givens rotations defining Q."
;;   (let ((m (nrows matrix))
;;         (n (ncols matrix))
;;         (gs nil))
;;     (assert (cl:= m n))
;;     (loop :for i :below (1- n)
;;           :for (c s) := (multiple-value-list
;;                          (givens-entries (tref matrix i i) (tref matrix (1+ i) i)))
;;           :for g := (make-givens-rotation :i i :j (1+ i) :c c :s s)
;;           :do (left-apply-givens! matrix g :start-idx i) ; TODO: end-idx = min(i+2, (1- n))
;;               (push g gs))
;;     gs))


;; (defun incf-diag! (matrix val)
;;   (assert (cl:= (nrows matrix) (ncols matrix)))
;;   (dotimes (i (nrows matrix))
;;     (incf (tref matrix i i) val))
;;   matrix)

(defun wilkinson-shift (A)
  (assert (cl:= (nrows A) (ncols A)))
  ;; bottom right 2x2 block looks like
  ;; [ ajj ajk ; akj akk]
  ;; The following is from Golub & Van Loan sec 8.3.5
  (let* ((k (1- (nrows A)))
         (j (1- k))
         (ajj (tref A j j))
         (ajk (tref A j k))
         (akk (tref A k k))
         (ajk2 (* ajk (conjugate ajk)))
         (delta (/ (- ajj akk) 2))
         (d2 (* delta delta))
         (sgnd (if (minusp delta) -1d0 1d0)) ; we do NOT want SIGNUM here
         )
    (- akk
       (/ ajk2
          (+ delta
             (* sgnd (sqrt (+ d2 ajk2))))))))

(defun qrstep! (A)
  (let* ((n (ncols A))
         (shift (wilkinson-shift A))
         (x (- (tref A 0 0) shift))
         (z (tref A 1 0)))
    (loop :for k :below (1- n)
          :for (c s) := (multiple-value-list (givens-entries x z))
          :for g := (make-givens-rotation :i k :j (1+ k) :c c :s s)
          :do (left-apply-givens! A g)
              (right-apply-givens! A g)
          :when (< k (- n 2))
            :do (setf x (tref A (1+ k) k)
                      z (tref A (+ 2 k) k))
          :collect g)))


(defun matrix-realpart (matrix)
  (typecase matrix
    (matrix/complex-single-float
     (let ((result (empty (shape matrix) :type 'single-float)))
      (map-to #'realpart matrix result)
       result))
    (matrix/complex-double-float
     (let ((result (empty (shape matrix) :type 'double-float)))
      (map-to #'realpart matrix result)
       result))
    (otherwise matrix)))


(defmethod hermitian-eig-lisp (matrix)
  (assert (hermitian-matrix-p matrix))
  (let ((m (ncols matrix))        
        (eigs nil)
        (matrix (deep-copy-tensor matrix)))
    (multiple-value-bind (H Q) (hessenberg matrix)
      (setf H (matrix-realpart H))
      (loop :with i := (1- m)
            :until (zerop i)
            :do (dolist (g (qrstep! H))
                  (right-apply-givens! Q g))
            :when (<= (abs (tref H (1- i) i))
                      (* *double-comparison-threshold*
                         (+ (abs (tref H i i)) (abs (tref H (1- i) (1- i))))))
              :do (push (tref H i i) eigs)
                  (setf H (slice H (list 0 0) (list i i)))
                  (decf i)
            :finally (push (tref H 0 0) eigs))
      (values eigs Q))))

(defun hilbert-matrix (n)
  (let ((H (empty (list n n) :type 'double-float)))
    (into! (lambda (i j) (/ (+ i j 1)))
           H)
    H))

(defun random-hermitian (n)
  (let ((a (rand (list n n) :type '(complex double-float)))
        (b (rand (list n n) :type '(complex double-float))))
    (scale! b #C(0d0 1d0))
    (let ((c (.+ a b)))
      (.+ c (conjugate-transpose c)))))
