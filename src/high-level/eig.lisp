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


(defstruct householder-reflection
  "Representation of a Householder reflection."
  v (idx 0))

(defun column->householder (A i0 j0)
  "Construct a Householder reflection which, when applied from the left, puts zeros in column J0 below row I0."
  (let* ((m (nrows A))      
         (v (zeros (list (- m i0)) :type (element-type A))))
    (loop :for i :from i0 :below m
          :for iv :from 0
          :do (setf (tref v iv) (tref A i j0)))
    (let ((norm-v (norm v)))
      (if (zerop norm-v)
          (make-householder-reflection)
          (let ((v0 (tref v 0)))
            (incf (tref v 0)                ; TODO: revisit this
                  (* norm-v (/ v0 (abs v0))))
            (scale! v (/ (norm v)))
            (make-householder-reflection :v v :idx j0))))))

(defun row->householder (A i0 j0)
  "Construct a Householder reflection which, when applied from the right, puts zeros in row I0 to the right of column J0."
  (let* ((n (ncols A))
         (v (zeros (list (- n j0)) :type (element-type A))))
    (loop :for j :from j0 :below n
          :for jv :from 0
          :do (setf (tref v jv) (conjugate (tref A i0 j))))
    (let ((norm-v (norm v)))
      (if (zerop norm-v)
          (make-householder-reflection)
          (let ((v0 (tref v 0)))
            (incf (tref v 0)
                  (* norm-v (/ v0 (abs v0))))
            (scale! v (/ (norm v)))
            (make-householder-reflection :v v :idx i0))))))

(defun left-apply-householder! (A hh)
  "Apply the Householder reflection HH to A[i:,j:], from the left."
  (when (null (householder-reflection-v hh))
    (return-from left-apply-householder! A))
  (let* ((m (nrows A))
         (n (ncols A))
         (v (householder-reflection-v hh))
         (i0 (- m (size v))))
    (flet ((column-reflect! (j)
             "Reflect A[i0:m,j]."
             (let ((v-dot-A
                     (loop :for i :from i0 :below m
                           :for iv :from 0
                           :sum (* (conjugate (tref v iv)) (tref A i j)))))
               (loop :for i :from i0 :below m
                     :for iv :from 0
                     :do (decf (tref A i j)
                               (* 2 (tref v iv) v-dot-A))))))
      (loop :for j :from (householder-reflection-idx hh) :below n
            :do (column-reflect! j))))
  A)

(defun right-apply-householder! (A hh)
  "Apply the Householder reflection HH to A[i:,j:] from the right."
  (declare (optimize debug))
  (when (null (householder-reflection-v hh))
    (return-from right-apply-householder! A))
  (let* ((m (nrows A))
         (n (ncols A))
         (v (householder-reflection-v hh))
         (j0 (- n (size v))))
    (flet ((row-reflect! (i)
             "Reflect A[i,j0:m] across (dagger V)."
             (let ((v*-dot-a
                     (loop :for j :from j0 :below n
                           :for iv :from 0
                           :sum (* (tref A i j) (tref v iv)))))
               (loop :for j :from j0 :below n
                     :for iv :from 0
                     :do (decf (tref A i j)
                               (* 2 v*-dot-a (conjugate (tref v iv))))))))
      (loop :for i :from (householder-reflection-idx hh) :below m
            :do (row-reflect! i))))
  A)

;;; TODO: - rename HESSENBERG to TRIDIAGONAL-FORM
;;;       - get Q matrix from above
;;;       - apply givens rotation to Q

(defun scale-column! (matrix idx scalar)
  (loop :for i :below (nrows matrix)
        :do (setf (tref matrix i idx)
                  (* scalar (tref matrix i idx)))
        :finally (return matrix)))

(defun scale-row! (matrix idx scalar)
  (loop :for j :below (ncols matrix)
        :do (setf (tref matrix idx j)
                  (* scalar (tref matrix idx j)))
        :finally (return matrix)))

;;; QR

(defun %qr-lisp (matrix &key (reduced t))
  (let* ((m (nrows matrix))
	 (n (ncols matrix))
         (p (min m n))
	 (hhs nil)
         (Q (eye (list m (if reduced p m)) :type (element-type matrix)))
	 (R (deep-copy-tensor matrix)))
    ;; compute R and vs
    (loop :for k :below p
          :for hh := (column->householder R k k)
          :do (left-apply-householder! R hh)
              (push hh hhs))
    ;; compute Q
    (loop :for hh :in hhs
          :do (left-apply-householder! Q hh))
    (when reduced
      ;; reduced factorization; trim R
      (setf R (slice R (list 0 0) (list p (ncols R)))))
    ;; TODO: fix this with the choice of householder coeffs
    ;; force positive values on the diagonal
    (loop :for i :below (min (nrows R) (ncols R))
          :for val := (tref R i i)
          :unless (and (zerop (imagpart val))
                       (not (minusp (realpart val))))
            :do (scale-column! Q i (/ val (abs val)))
                (scale-row! R i (/ (conjugate val) (abs val))))
    (values Q R)))

(defmethod qr-lisp ((matrix matrix))
  (%qr-lisp matrix :reduced (< (ncols matrix) (nrows matrix))))

;;;

(defun hessenberg (matrix)
  "Reduce MATRIX to Hessenberg form via a sequence of unitary similarity transformations."
  (let* ((m (nrows matrix))
	 (n (ncols matrix))
	 (hhs nil)
	 (H (deep-copy-tensor matrix))
         (Q (eye (list m n) :type (element-type matrix))))
    (loop :for k :below (1- m)
          :for hh := (column->householder H (1+ k) k)
          :do (left-apply-householder! H hh)
              (right-apply-householder! H hh)
              (push hh hhs))
    ;; update Q
    (loop :for hh :in hhs
          :do (left-apply-householder! Q hh))
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
    (values H Q)))


(defun bidiagonal (matrix)
  "Reduce MATRIX to bidiagonal form via left and right multiplication by unitary transformations."
  (declare (optimize debug))
  (let* ((m (nrows matrix))
	 (n (ncols matrix))
         (left-hhs nil)
	 (right-hhs nil)
	 (H (deep-copy-tensor matrix)))
    (assert (cl:= m n))
    (dotimes (k m)
      (let ((hh (column->householder H k k)))
        (left-apply-householder! H hh)
        (push hh left-hhs))
      (when (< k (- n 2))
        (let ((hh (row->householder H k (1+ k))))
          (right-apply-householder! H hh)
          (push hh right-hhs))))
    (let ((U (eye (list m m) :type (element-type matrix)))
          (V (eye (list n n) :type (element-type matrix))))
      (loop :for hh :in left-hhs
            :do (left-apply-householder! U hh))
      (loop :for hh :in right-hhs
            :do (left-apply-householder! V hh))
      (values U H V))))


(defun wilkinson-shift (ajj ajk akk)
  "Get the eigenvalue of [ajj ajk; ajk akk] closest to akk."
  ;; The following is from Golub & Van Loan sec 8.3.5
  (let* ((ajk2 (* ajk (conjugate ajk)))
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
         (shift (wilkinson-shift (tref A (- n 2) (- n 2))
                                 (tref A (- n 2) (1- n))
                                 (tref A (1- n) (1- n))))
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

(defun svdstep! (B)
  (flet ((abs2 (x)
           (* (conjugate x) x)))
    (let* ((m (nrows B))
           (n (ncols B))
           (Bij (if (> m 2) (tref B (- m 3) (- n 2)) 0))
           (Bjj (tref B (- m 2) (- n 2)))
           (Bjk (tref B (- m 2) (1- n)))
           (Bkk (tref B (1- n) (1- n)))
           ;; we want to shift by the amount associated with the tridiagonal B^H B
           ;; the following can be worked out by hand or found in 8.6 of Golub & van Loan
           (shift (wilkinson-shift (+ (abs2 Bij) (abs2 Bjj))
                                   (* Bjj Bjk)
                                   (+ (abs2 Bjk) (abs2 Bkk))))
           (y (- (tref B 0 0) shift))
           (z (tref B 0 1))
           (gs nil))
      (dotimes (k (1- n))
        (multiple-value-bind (c s) (givens-entries y z)
          (let ((g (make-givens-rotation :i k :j (1+ k) :c c :s s)))
            (right-apply-givens! B g)
            (push g gs)
            (setf y (tref B k k)
                  z (tref B (1+ k) k)))
          (multiple-value-bind (c s) (givens-entries y z)
            (let ((g (make-givens-rotation :i k :j (1+ k) :c c :s s)))
              (left-apply-givens! B g)
              (push g gs)
              (when (< k (- n 2))
                (setf y (tref B k (1+ k))
                      z (tref B k (+ k 2))))))))
      (nreverse gs))))

;;; TODO: consistent use of M and N


(defun unitary-extension (columns)
  (let ((m (nrows columns))
        (n (ncols columns)))
    (assert (<= n m))
    (let ((u (zeros (list m m) :type (element-type columns))))
      (slice-to columns (list 0 0) (list m n)
                u (list 0 0))
      (multiple-value-bind (q r) (qr u)
        (declare (ignore r))
        q))))

(defun svd-lisp (matrix &key reduced)
  ;; short and fat => find svd of tranpose
  (when (> (ncols matrix) (nrows matrix))
    (multiple-value-bind (U D Vh) (svd-lisp (transpose matrix) :reduced reduced)
      (return-from svd-lisp
        (values (transpose! Vh :fast t)
                (transpose! D :fast t)
                (transpose! U :fast t)))))
  (multiple-value-bind (Q R) (%qr-lisp matrix)
    (let ((svals nil))
      (multiple-value-bind (U B V) (bidiagonal R)
        ;; TODO: realpart
        (assert (cl:= (nrows b) (ncols b)))
        (loop :with i := (1- (nrows B))
              :until (zerop i)
              :do (loop :for (gv gu) :on (svdstep! B) :by #'cddr
                        :do (right-apply-givens! U gu)
                            (right-apply-givens! V gv))
              :when (<= (abs (tref B (1- i) i))
                        (* 2 *double-comparison-threshold* ; TODO
                           (+ (abs (tref B (1- i) (1- i))) (abs (tref B i i)))))
                :do (push (tref B i i) svals)
                    (setf B (slice B (list 0 0) (list i i)))
                    (decf i)
              :finally (push (tref B 0 0) svals))
        (setf U (@ Q U))
        ;; permute singular values
        ;; we want positive values, from greatest to smallest
        (let ((U-sorted (zeros (shape U) :type (element-type U)))
              (Vt-sorted (zeros (shape V) :type (element-type V))))
          (let ((svals-sorted
                  (loop :with tagged := (loop :for i :from 0 :for v :in svals :collect (cons i v))
                        :for (old-idx . val) :in (sort tagged #'> :key (lambda (x) (abs (cdr x))))
                        :for new-idx :from 0
                        :do (dotimes (i (nrows U))
                              (setf (tref U-sorted i new-idx)                                    
                                    (if (minusp val)
                                        (- (tref U i old-idx))
                                        (tref U i old-idx))))
                            (dotimes (j (ncols V))
                              (setf (tref Vt-sorted j new-idx)
                                    (tref V old-idx j)))
                        :collect (abs val))))
            (if reduced
                (values U-sorted
                        (from-diag svals-sorted :type (element-type matrix))
                        Vt-sorted)
                (let ((d (zeros (list (nrows matrix) (ncols matrix)) :type (element-type matrix))))
                  (loop :with etype := (element-type matrix)
                        :for i :below (min (nrows matrix) (ncols matrix))
                        :do (setf (tref d i i) (coerce (pop svals-sorted) etype)))
                  (values (unitary-extension U-sorted)
                          d
                          Vt-sorted)))))))))


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
