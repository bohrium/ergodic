module Main where

import Graphics.Gloss
import Graphics.Gloss.Interface.Environment
import System.Random

-- x, y, dx, dy
type Particle
        = (Float, Float, Float, Float)



nb::Int
nb = 128
half::Int
half = nb `div` 2

youngs::Float
youngs = 0.0001
lambda::Float
lambda = 0.5


main :: IO ()
main
 = do   g <- getStdGen
        (width,height) <- getScreenSize
        let initialstate = generateParticles g width height
        simulate window background fps initialstate render update
 where
        window          = FullScreen
        background      = black
        fps             = 30
        render xs       = (pictures $
                            (map (particleImage white) (take nb xs)) ++
                            (map (particleImage blue) (take half $ drop nb xs))++ 
                            (map (particleImage red) (drop (nb+half) xs))) 
        update _ dt     = updateParticles 0.20

toF :: Int -> Float
toF n = fromIntegral n::Float

-- | Generates particles from StdGen
generateParticles :: StdGen -> Int -> Int -> [Particle]
generateParticles gen widthInt heightInt
 = map (g . f) [0 .. (nb-1+half+half)]
 where
        f = \n -> ((fromIntegral n::Float) / (fromIntegral nb::Float) - 0.5,
                   0.0, 0.0, 0.00*cos ((toF n) * 2.0 * (pi::Float)/(toF nb))
                             + (if 5<=n && n<10 then 0.01 else 0.0) - (10-5)*0.01/(toF nb)  
                             + (if 8<=n && n<9  then 0.01 else 0.0) - ( 9-8)*0.01/(toF nb)) 
        g = \(x,y,p_x,p_y) -> (
                x*1,y*1,p_x*1,p_y*1
            )

get_y :: Particle -> Float
get_y (_,y,_,_) = y
get_p :: Particle -> Float
get_p (_,_,_,p) = p

fourier :: Int -> [Particle] -> Particle
fourier k ps
    = (0.01*kF, -0.5+2*sqrt((yya^2+yyb^2)*kF*kF+(ppa^2+ppb^2)/youngs), 0, 0)
    where
        yya = (sum $ map (\n -> (cos ((toF k) * (toF n) * scale))*(get_y (ps!!n))) [0 .. (nb-1)]) / (toF nb)
        yyb = (sum $ map (\n -> (sin ((toF k) * (toF n) * scale))*(get_y (ps!!n))) [0 .. (nb-1)]) / (toF nb)
        ppa = (sum $ map (\n -> (cos ((toF k) * (toF n) * scale))*(get_p (ps!!n))) [0 .. (nb-1)]) / (toF nb)
        ppb = (sum $ map (\n -> (sin ((toF k) * (toF n) * scale))*(get_p (ps!!n))) [0 .. (nb-1)]) / (toF nb)
        scale = 2*(pi::Float)/(toF nb)
        kF = toF k

-- | Particle to its picture
particleImage :: Color -> Particle -> Picture
particleImage c (x,y,_,_)
 = translate (960*x) (960*y) $ color c $ circleSolid 2

-- | To update particles for next frame
updateParticles :: Float -> [Particle] -> [Particle]
updateParticles half_dt ps
 = physicals ++ now_freqs ++ avg_freqs
    where physicals = (update_atom . update_atom . update_atom) (take nb ps)
          update_atom = (  (accelerateParticles    half_dt) 
                         . (moveParticles       (2*half_dt)) 
                         . (accelerateParticles    half_dt)) 
          now_freqs = (map (\k -> fourier k physicals) [0 .. (half-1)])
          avg_freqs = [(x, y+0.01*((y_new+0.25)-y),0,0) | ((x,y_new,_,_),(_,y,_,_)) <- take half $ zip (drop nb ps) (drop (nb + half) ps)]

-- | Moves particles based on their speed
moveParticles :: Float -> [Particle] -> [Particle]
moveParticles dt
 = map (\(x,y,dx,dy) -> (x+dx*dt,y+dy*dt,dx,dy))

-- | Accelerates particles based on gravity
accelerateParticles :: Float -> [Particle] -> [Particle]
accelerateParticles dt ps
 = [accelerated_at dt ps i    | i <- [0..(nb-1)]]

accelerated_at :: Float -> [Particle] -> Int -> Particle
accelerated_at dt ps i
    = (x,ya,px,py+dt*(force_from yz ya yb))
      where (_,yz, _, _) = ps!!((i-1) `mod` nb)
            (x,ya,px,py) = ps!!( i    `mod` nb)
            (_,yb, _, _) = ps!!((i+1) `mod` nb)

force_from :: Float -> Float -> Float -> Float 
force_from z a b =
    (youngs*(toF nb)*(pi::Float)) * 
            ((z-a) * (1.0 + lambda * ((z-a)*(toF nb))^2) + 
             (b-a) * (1.0 + lambda * ((b-a)*(toF nb))^2)) 
