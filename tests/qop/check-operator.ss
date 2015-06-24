(define (list-equal? la lb)
  (cond
   [(null? la) (null? lb)]
   [(null? lb) #f]
   [(= (car la) (car lb)) (list-equal? (cdr la) (cdr lb))]
   [else #f]))

(define (x-ratio a b)
  (if (< (abs b) 1e-10) 0.0 (exact->inexact (/ a b))))

(define (get-scale f1 f2)
  (for-each (lambda (a b)
              (if (not (list-equal? (car a) (car b)))
                  (printf "************ site mismatch: ~a ~a~%" (car a) (car b))
                  (printf "~20,7e  ~20,7e -- ~a~%"
                          (x-ratio (cadr a) (cadr b))
                          (x-ratio (caddr a) (caddr b))
                          (car a))))
            f1 f2))
