module Main where

import Graphics.Gloss
import Graphics.Gloss.Interface.Environment
import System.Random

-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-- ~~~~~~~~  0. Simulation Parameters  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-- (d/dt)^2 x = Y 2 (cos(tau k/N)-1) x
-- x = C exp(tau t/T) --> tau T = sqrt(2Y (1 - cos(tau k/N)))
-- energy = (N/2) * 2Y (cos(tau k/N)-1)^2 * amplitude^2/2

enum xs = zip [0..] xs

type Particle = (Float, Float) -- Position q, Momentum p
get_q :: Particle -> Float
get_p :: Particle -> Float
get_q (q,_) = q
get_p (_,p) = p

type RenderableState = (Int, [Particle], [Float], [Float])
get_qs     ::RenderableState->[Float]
get_fourier::RenderableState->[Float]
get_avg    ::RenderableState->[Float]
get_i       (i,qps,_,_) = i
get_qs      (i,qps,_,_) = map get_q qps
get_fourier (i,_  ,f,_) = f
get_avg     (i,_  ,_,a) = a

nb::Int
nb = 64

half::Int
half = nb `div` 2

tau::Float
tau = 2*pi
costau::Float->Float
sintau::Float->Float
costau z = cos (tau*z) 
sintau z = sin (tau*z) 

konst::Float
konst = 0.5 
lambda::Float
lambda = 1.0

time_step::Float
time_step = 0.10

-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-- ~~~~~~~~  1. Rendering and Main Loop  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

main :: IO ()
main
  = do g <- getStdGen
       (width, height) <- getScreenSize
       let initial_state = (0,physicals,
                            replicate half 0.0,
                            getFour physicals)
            where physicals = generateParticles g width height
       simulate window background fps initial_state render update
    where window      = FullScreen
          background  = black
          fps         = 60
          render state = (pictures $
                          (graphImage white  0.0 0.05 (-0.75) 0.25 $ get_qs      state) ++
                          (graphImage violet 0.0 0.05 ( 0.50) 4.00 $ replicate half 0.0)++
                          (graphImage blue   0.0 0.05 ( 0.50) 4.00 $ map sqrt $ get_fourier state) ++ 
                          (graphImage violet 0.0 0.05 ( 1.20) 4.00 $ replicate half 0.0)++
                          (graphImage red    0.0 0.05 ( 1.20) 4.00 $ map sqrt $ get_avg     state))
          update _ _ = updateParticles time_step

toF :: Int -> Float
toF n = fromIntegral n::Float

-- | initial conditions
generateParticles :: StdGen -> Int -> Int -> [Particle]
generateParticles gen widthInt heightInt
  = map f [0 .. (nb-1)]
    where f n = ((sum $ map (\k -> (1.5/(toF k)) * (cos_scale k n)           ) [1..10]),  
                 (sum $ map (\k -> (1.5/(toF k)) * (sin_scale k n) * (pamp k)) [1..10])) 
          cos_scale k n = costau ((toF k)*(toF n)/nbF) 
          sin_scale k n = sintau ((toF k)*(toF n)/nbF) 
          nbF = toF nb
          pamp k = sqrt $ konst * ((sin_scale k 1)^2 + (1 - cos_scale k 1)^2)

fourier :: Int -> [Particle] -> Float
fourier k qps
  = eff_konstant*(qqa^2+qqb^2) + (ppa^2+ppb^2)
    where
        qqa = (sum $ map (\(n,(q,p)) -> (cos_scale n) * q) $ enum qps) / nbF
        qqb = (sum $ map (\(n,(q,p)) -> (sin_scale n) * q) $ enum qps) / nbF
        ppa = (sum $ map (\(n,(q,p)) -> (cos_scale n) * p) $ enum qps) / nbF
        ppb = (sum $ map (\(n,(q,p)) -> (sin_scale n) * p) $ enum qps) / nbF
        cos_scale n = costau ((toF k)*(toF n)/nbF) 
        sin_scale n = sintau ((toF k)*(toF n)/nbF) 
        nbF = toF nb 
        eff_konstant = konst * ((sin_scale 1)^2 + (1 - cos_scale 1)^2)

-- | graph to its pictures
graphImage :: Color -> Float -> Float -> Float -> Float -> [Float] -> [Picture]
graphImage c x dx y dy fs
 = map (\(n,f) -> translate (pix*(horiz n)) (pix*(verti f)) dot) $ enum fs
   where horiz n = x + dx * ((toF n)-(toF $ length fs)/2) 
         verti f = y + dy * f 
         dot = color c $ circleSolid 3
         pix = 240 

-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-- ~~~~~~~~  2. Physical Dynamics  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-- | to update particles for next frame
updateParticles :: Float -> RenderableState -> RenderableState
updateParticles half_dt (i,qps,fs,as)
 = (i+1,physicals,now_freqs,avg_freqs) 
    where physicals = super_update qps
          update_atom = (  (accelerateParticles    half_dt) 
                         . (moveParticles       (2*half_dt)) 
                         . (accelerateParticles    half_dt)) 
          now_freqs = getFour physicals 
          avg_freqs = [a + 0.002 * (f-a) | (f,a) <- zip fs as]
          super_update = if i `mod` 2 == 0
                         then (update_atom . update_atom . update_atom . update_atom . update_atom) 
                         else (update_atom . update_atom . update_atom . update_atom) 

getFour physicals = (map (\k -> fourier k physicals) [0 .. (half-1)])

-- | moves particles based on their speed
moveParticles :: Float -> [Particle] -> [Particle]
moveParticles dt
 = map (\(q, p) -> (q + p*dt, p))

-- | accelerates particles based on spring forces
accelerateParticles :: Float -> [Particle] -> [Particle]
accelerateParticles dt ps
 = [accelerated_at dt ps i    | i <- [0..(nb-1)]]

accelerated_at :: Float -> [Particle] -> Int -> Particle
accelerated_at dt qps i
    = (qa, pa+(force_from qz qa qb)*dt)
      where (qz, _) = qps!!((i-1) `mod` nb)
            (qa,pa) = qps!!( i    `mod` nb)
            (qb, _) = qps!!((i+1) `mod` nb)

force_from :: Float -> Float -> Float -> Float 
force_from z a b  
  = konst * (za * (1.0 + lambda * za^2) + 
             ba * (1.0 + lambda * ba^2)) 
    where za = (z-a)
          ba = (b-a) 
