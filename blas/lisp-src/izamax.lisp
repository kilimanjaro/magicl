;;; Compiled by f2cl version:
;;; ("f2cl1.l,v 95098eb54f13 2013/04/01 00:45:16 toy $"
;;;  "f2cl2.l,v 95098eb54f13 2013/04/01 00:45:16 toy $"
;;;  "f2cl3.l,v 96616d88fb7e 2008/02/22 22:19:34 rtoy $"
;;;  "f2cl4.l,v 96616d88fb7e 2008/02/22 22:19:34 rtoy $"
;;;  "f2cl5.l,v 95098eb54f13 2013/04/01 00:45:16 toy $"
;;;  "f2cl6.l,v 1d5cbacbb977 2008/08/24 00:56:27 rtoy $"
;;;  "macros.l,v 1409c1352feb 2013/03/24 20:44:50 toy $")

;;; Using Lisp SBCL 2.0.9
;;; 
;;; Options: ((:prune-labels nil) (:auto-save t) (:relaxed-array-decls t)
;;;           (:coerce-assigns :as-needed) (:array-type ':array)
;;;           (:array-slicing t) (:declare-common nil)
;;;           (:float-format single-float))

(in-package :magicl.lisp-blas)


(defun izamax (n zx incx)
  (declare (type (array f2cl-lib:complex16 (*)) zx)
           (type (f2cl-lib:integer4) incx n))
  (f2cl-lib:with-multi-array-data
      ((zx f2cl-lib:complex16 zx-%data% zx-%offset%))
    (prog ((i 0) (ix 0) (dmax 0.0d0) (izamax 0))
      (declare (type (double-float) dmax)
               (type (f2cl-lib:integer4) izamax ix i))
      (setf izamax 0)
      (if (or (< n 1) (<= incx 0))
          (go end_label))
      (setf izamax 1)
      (if (= n 1)
          (go end_label))
      (cond
       ((= incx 1)
        (setf dmax (dcabs1 (f2cl-lib:fref zx-%data% (1) ((1 *)) zx-%offset%)))
        (f2cl-lib:fdo (i 2 (f2cl-lib:int-add i 1))
                      ((> i n) nil)
          (tagbody
            (cond
             ((> (dcabs1 (f2cl-lib:fref zx (i) ((1 *)))) dmax) (setf izamax i)
              (setf dmax
                      (dcabs1
                       (f2cl-lib:fref zx-%data% (i) ((1 *)) zx-%offset%)))))
           label100000)))
       (t (setf ix 1)
        (setf dmax (dcabs1 (f2cl-lib:fref zx-%data% (1) ((1 *)) zx-%offset%)))
        (setf ix (f2cl-lib:int-add ix incx))
        (f2cl-lib:fdo (i 2 (f2cl-lib:int-add i 1))
                      ((> i n) nil)
          (tagbody
            (cond
             ((> (dcabs1 (f2cl-lib:fref zx (ix) ((1 *)))) dmax) (setf izamax i)
              (setf dmax
                      (dcabs1
                       (f2cl-lib:fref zx-%data% (ix) ((1 *)) zx-%offset%)))))
            (setf ix (f2cl-lib:int-add ix incx))
           label100001))))
      (go end_label)
     end_label
      (return (values izamax nil nil nil)))))

(in-package #-gcl #:cl-user #+gcl "CL-USER")
#+#.(cl:if (cl:find-package '#:f2cl) '(and) '(or))
(eval-when (:load-toplevel :compile-toplevel :execute)
  (setf (gethash 'fortran-to-lisp::izamax
                 fortran-to-lisp::*f2cl-function-info*)
          (fortran-to-lisp::make-f2cl-finfo :arg-types
                                            '((fortran-to-lisp::integer4)
                                              (array fortran-to-lisp::complex16
                                               (*))
                                              (fortran-to-lisp::integer4))
                                            :return-values '(nil nil nil)
                                            :calls '(fortran-to-lisp::dcabs1))))
