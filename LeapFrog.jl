#This library contains the "cross" function used for cross product calculations
using LinearAlgebra
#This library containts the necessary functions for graphing, specifically scatter plots
using Plots

#=
The leapFrog function outputs a plot showing the varrying x and y coordinates of
four vortecies that start equadistant from each other in the shape of a square. The
top two verticies rotate counterclockwise while the bottom two rotate clockwise.
The inputs necessary for the function to work are as follows
1) Circulation in the form of an 1X3 array, [x,y,z]
2) Initial distance betweeb vortices,i.e. side length of the square.
3) Time step
4) Total time
=#

function LeapFrog(Gamma,d,dt,t)

#Calculates the number of data points for the time steps over the time period
itterations = Int(t/dt + 1);

#Creates arrays for for the velocities and positions of the four vortices
V1 = zeros(itterations,3);
V2 = zeros(itterations,3);
V3 = zeros(itterations,3);
V4 = zeros(itterations,3);

P1 = zeros(itterations,3);
P2 = zeros(itterations,3);
P3 = zeros(itterations,3);
P4 = zeros(itterations,3);

#Creates arrays to store the distances between all the vortices for each itteration
P1toP2 = zeros(itterations,3);
P1toP3 = zeros(itterations,3);
P1toP4 = zeros(itterations,3);

P2toP1 = zeros(itterations,3);
P2toP3 = zeros(itterations,3);
P2toP4 = zeros(itterations,3);

P3toP1 = zeros(itterations,3);
P3toP2 = zeros(itterations,3);
P3toP4 = zeros(itterations,3);

P4toP1 = zeros(itterations,3);
P4toP2 = zeros(itterations,3);
P4toP3 = zeros(itterations,3);

#Sets the initial position of each vortex symetrically about the x axis, starting at x=0.
P1[1,:] = [0,-d/2,0];
P2[1,:] = [0,d/2,0];
P3[1,:] = [d,d/2,0];
P4[1,:] = [d,-d/2,0];

#Creates counter to enter while loop. Starts at 2 since initial positions were previously set, and velocity starts at zero
i = 2;

#While loop to calculate velocity and position for each vortex in each itteration
while i <= itterations

    #Calculates the distances between all the vortices
    P1toP2[i-1,:] = P2[i-1,:]-P1[i-1,:];
    P1toP3[i-1,:] = P3[i-1,:]-P1[i-1,:];
    P1toP4[i-1,:] = P4[i-1,:]-P1[i-1,:];

    P2toP1[i-1,:] = -1*P1toP2[i-1,:];
    P2toP3[i-1,:] = P3[i-1,:]-P2[i-1,:];
    P2toP4[i-1,:] = P4[i-1,:]-P2[i-1,:];

    P3toP1[i-1,:] = -1*P1toP3[i-1,:];
    P3toP2[i-1,:] = -1*P2toP3[i-1,:];
    P3toP4[i-1,:] = P4[i-1,:]-P3[i-1,:];

    P4toP1[i-1,:] = -1*P1toP4[i-1,:];
    P4toP2[i-1,:] = -1*P2toP4[i-1,:];
    P4toP3[i-1,:] = -1*P3toP4[i-1,:];

    #Sums the relative velocities between vorticies. Note that counterclockwise vorticities cause a negative velocity
    V1[i,:] = cross(Gamma,P1toP4[i-1,:])/(2*pi*(norm(P1toP4[i-1,:]))^2) - cross(Gamma,P1toP2[i-1,:])/(2*pi*(norm(P1toP2[i-1,:]))^2) - cross(Gamma,P1toP3[i-1,:])/(2*pi*(norm(P1toP3[i-1,:]))^2);

    V2[i,:] = cross(Gamma,P2toP1[i-1,:])/(2*pi*(norm(P2toP1[i-1,:]))^2) + cross(Gamma,P2toP4[i-1,:])/(2*pi*(norm(P2toP4[i-1,:]))^2) - cross(Gamma,P2toP3[i-1,:])/(2*pi*(norm(P2toP3[i-1,:]))^2);

    V3[i,:] = cross(Gamma,P3toP1[i-1,:])/(2*pi*(norm(P3toP1[i-1,:]))^2) + cross(Gamma,P3toP4[i-1,:])/(2*pi*(norm(P3toP4[i-1,:]))^2) - cross(Gamma,P3toP2[i-1,:])/(2*pi*(norm(P3toP2[i-1,:]))^2);

    V4[i,:] = cross(Gamma,P4toP1[i-1,:])/(2*pi*(norm(P4toP1[i-1,:]))^2) - cross(Gamma,P4toP2[i-1,:])/(2*pi*(norm(P4toP2[i-1,:]))^2) - cross(Gamma,P4toP3[i-1,:])/(2*pi*(norm(P4toP3[i-1,:]))^2);

    #Calculates the new position by using xf = xo + Vdt
    P1[i,:] = P1[i-1,:] + dt*V1[i,:];

    P2[i,:] = P2[i-1,:] + dt*V2[i,:];

    P3[i,:] = P3[i-1,:] + dt*V3[i,:];

    P4[i,:] = P4[i-1,:] + dt*V4[i,:];

    #Rounds position coordinates to the nearest thousandth place for the the first few calculations
    if i<5

       P1[i,1] = round(P1[i,1], digits = 3);
       P1[i,2] = round(P1[i,2], digits = 3);
       P1[i,3] = round(P1[i,3], digits = 3);

       P2[i,1] = round(P2[i,1], digits = 3);
       P2[i,2] = round(P2[i,2], digits = 3);
       P2[i,3] = round(P2[i,3], digits = 3);

       P3[i,1] = round(P3[i,1], digits = 3);
       P3[i,2] = round(P3[i,2], digits = 3);
       P3[i,3] = round(P3[i,3], digits = 3);

       P4[i,1] = round(P4[i,1], digits = 3);
       P4[i,2] = round(P4[i,2], digits = 3);
       P4[i,3] = round(P4[i,3], digits = 3);

    end

    i = i+1;

end

#Plots the position arrays
plot(P1[:,1], P1[:,2], seriestype = :scatter, title = "Leap Frog", markerstrokecolor = :blue,markercolor = :blue)
plot!(P2[:,1], P2[:,2], seriestype = :scatter, markerstrokecolor = :blue,markercolor = :blue)
plot!(P3[:,1], P3[:,2], seriestype = :scatter, markerstrokecolor = :orange, markercolor = :orange)
plot!(P4[:,1], P4[:,2], seriestype = :scatter, markerstrokecolor = :orange, markercolor = :orange, legend = false, xaxis = ("Distance"), yaxis = ("Radial Location"))

savefig(leepfrog)

end
