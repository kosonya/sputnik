#!python
# -*- coding:utf8 -*-
# Copyright (C) 2024 Sophia Kovaleva
# sophia.m.kovaleva@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pygame
import sys
import random
import math
import numpy as np

G_constant = 6.674e-11 #N*m2*kg-2

def vector_distance(vec1, vec2):
    return math.sqrt(sum([(vec1[_x] - vec2[_x])**2 for _x in range(min(len(vec1), len(vec2)))]))


class Background:
    def __init__(self, screen_size, star_count=5000, star_color=(200, 200, 200)):
        #Make it internally bigger for easier resizing
        self.screen_size = [_x*5 for _x in screen_size]
        
        self.surface = pygame.Surface(self.screen_size)
        self.surface.fill((0,0,0))
        for _ in range(star_count):
            center = random.randrange(1, self.screen_size[0]), random.randrange(1, self.screen_size[1])
            pygame.draw.circle(self.surface, star_color, center, 1)
    
    def render(self, surface, dest=(0,0), pixels_per_km=1):
        surface.blit(self.surface, dest)
    
    def get_size(self):
        return self.surface.get_size()


class Planet:
    def __init__(self, radius_km=6.37e3, mass_kg=5.97e24, atmosphere_altitude_km = 1000, default_pixels_per_km=1.0/32.0, fill_color=(255,255,255), ocean_color=(0, 0, 255), spot_color=(0, 255, 0), line_color=(0,0,0), n_spots=70, max_spot_size=0.25):
        self.radius_km = radius_km
        self.mass_kg = mass_kg
        self.center_coordinates_km = [0,0]
        self.pixels_per_km = default_pixels_per_km
        self.fill_color = fill_color
        self.ocean_color = ocean_color
        self.scale = 1
        self.spot_color = spot_color
        self.n_spots = n_spots
        self.max_spot_size = max_spot_size
        self.line_color = line_color
        self.atmosphere_altitude_km = atmosphere_altitude_km
        self.make_image()
        
        
    def make_image(self):
        self.planet_r_pixels = int(self.radius_km * self.pixels_per_km)
        self.atmosphere_pixels = int(self.atmosphere_altitude_km * self.pixels_per_km)
        self.image_size = [(self.planet_r_pixels+self.atmosphere_pixels)*2+2 for _ in [0,1]]
        self.image_center = [int(self.image_size[_x]/2) for _x in [0,1]]
        
        #Background
        self.surface = pygame.Surface(self.image_size)
        self.surface.fill(self.fill_color)
        self.surface.set_colorkey(self.fill_color)
        
        #Planet disk
        pygame.draw.circle(self.surface, self.ocean_color, self.image_center, self.planet_r_pixels)
        
        #Planet spots
        for _ in range(self.n_spots):
            r_spot = random.randrange(2, int(self.planet_r_pixels * self.max_spot_size))
            dr = random.randrange(1, int(self.planet_r_pixels - r_spot))
            vect = random.random()*2*math.pi
            x = int(self.image_center[0] + dr*math.cos(vect))
            y = int(self.image_center[1] + dr*math.sin(vect))
            pygame.draw.circle(self.surface, self.spot_color, (x, y), r_spot)
            
        #Planet grid
        pygame.draw.line(self.surface, self.line_color, (self.image_center[0]-self.planet_r_pixels, self.image_center[1]), (self.image_center[0]+self.planet_r_pixels, self.image_center[1]), 2)
        planet_rect = pygame.Rect(self.image_center[0]-self.planet_r_pixels, self.image_center[1]-self.planet_r_pixels, 2*self.planet_r_pixels, 2*self.planet_r_pixels)
        for w in [4, 2, 2, 2]:
            planet_rect.width -= self.planet_r_pixels/w
            planet_rect.centerx = self.image_center[0]
            pygame.draw.ellipse(self.surface, self.line_color, planet_rect, 2)
            
        #Planet atmosphere
        for altitude_pixels in range(self.atmosphere_pixels):
            color_gb = int(255.0/(altitude_pixels/3 + 1))
            pygame.draw.circle(self.surface, (0,color_gb,color_gb), self.image_center, self.planet_r_pixels + altitude_pixels, 1)
            
        
    def render(self, surface, dest=(0,0), pixels_per_km = 1):
        self.scale = pixels_per_km / self.pixels_per_km
        new_size = [int(_s * self.scale) for _s in self.surface.get_size()]
        out_surface = pygame.transform.scale(self.surface, new_size)
        out_surface.set_colorkey(self.fill_color)
        surface.blit(out_surface, dest)
        
    def get_size(self):
        return [int(_s * self.scale) for _s in self.surface.get_size()]
        
    def gravity_acceleration_ms2(self, location_km):
        global G_constant
        cen_cor_m = [_x*1000 for _x in self.center_coordinates_km]
        loc_m = [_x*1000 for _x in location_km]
        dist_m = vector_distance(cen_cor_m, loc_m)
        return G_constant * self.mass_kg / (dist_m**2)
        
    def gravity_acceleration_components_ms2(self, location_km):
        dx_km, dy_km = [c - l for (c, l) in zip(self.center_coordinates_km, location_km)]
        dist_km = vector_distance(self.center_coordinates_km, location_km)
        acc_ms2 = self.gravity_acceleration_ms2(location_km)
        acc_x_ms2 = acc_ms2 * dx_km / dist_km
        acc_y_ms2 = acc_ms2 * dy_km / dist_km
        return acc_x_ms2, acc_y_ms2

class Sputnik:
    def __init__(self, orbiting_planet, start_altitude_km = 1000, start_velocity_ms = 8000, mass_kg = 1000, thrust_kn = 100, default_pixels_per_km=1.0/32.0, color=(255, 0, 255), fill_color=(255,255,255), screen_size=[100, 100]):
        self.planet = orbiting_planet
        self.location_km = [0.0, start_altitude_km + self.planet.radius_km]
        self.velocity_ms = [start_velocity_ms, 0.0]
        self.acceleration_ms2 = [0.0, 0.0]
        self.mass_kg = mass_kg
        self.thrust_kn = thrust_kn
        self.color = color
        self.fill_color = fill_color
        self.pixels_per_km = default_pixels_per_km
        self.scale = 1
        self.radius_pixels = 10
        self.screen_size = screen_size
        self.image_size = self.screen_size
        self.image_center = [int(self.image_size[_x]/2) for _x in [0,1]]
        self.future = []
        self.future_color = (127, 127, 0)
        self.engine_firing = [0, 0]
        
        self.surface = pygame.Surface(self.image_size)
        self.surface.fill(self.fill_color)
        self.surface.set_colorkey(self.fill_color)
        
    def coordinates_to_pixels(self, location_km, pixels_per_km):
        #Assume that the planet is always at the center of the screen
        #Note that the Y axis is flipped for pixels
        pixel_coordinates = [
                    self.image_center[0] + int((location_km[0]+self.planet.center_coordinates_km[0])*pixels_per_km),
                    self.image_center[1] - int((location_km[1]+self.planet.center_coordinates_km[1])*pixels_per_km)
        ]
        return pixel_coordinates
        
    def render(self, surface, dest=(0,0), pixels_per_km=1):
        self.surface = pygame.Surface(self.image_size)
        self.surface.fill(self.fill_color)
        self.surface.set_colorkey(self.fill_color)
        
        #Assume that the planet is always at the center of the screen
        #Note that the Y axis is flipped for pixels
        self.pixels_per_km = pixels_per_km
        self.image_size = self.screen_size
        self.image_center = [int(self.image_size[_x]/2) for _x in [0,1]]
        self.pixel_coordinates = self.coordinates_to_pixels(self.location_km, pixels_per_km)

        if len(self.future) >= 2:
            pygame.draw.lines(self.surface, self.future_color, False, [self.coordinates_to_pixels(point, pixels_per_km) for (point, _) in self.future])
        
        pygame.draw.line(self.surface, (255,0,0), self.pixel_coordinates, [c-e*self.radius_pixels*2 for (c, e) in zip(self.pixel_coordinates, self.engine_firing)], 3)
        
        pygame.draw.circle(self.surface, self.color, self.pixel_coordinates, self.radius_pixels)
        
        surface.blit(self.surface, dest)

        
    def get_size(self):
        return self.image_size
        
    def compute_simulation_step(self, location_km, velocity_ms, acceleration_ms2, timestep_s):
        location_m = [_x*1000 for _x in location_km]
        new_velocity_ms = [v+a*timestep_s for (v, a) in zip(velocity_ms, acceleration_ms2)]
        new_location_km = [(s+v*timestep_s)/1000 for (s, v) in zip(location_m, velocity_ms)]
        return new_location_km, new_velocity_ms
        
    def make_simulation_step(self, timestep_s = 1.0):
        if self.engine_firing[0] != 0 or self.engine_firing[1] != 0:
            self.future = []
            self.acceleration_ms2 = self.planet.gravity_acceleration_components_ms2(self.location_km)
            self.acceleration_ms2 = [a+(e*self.thrust_kn*1000/self.mass_kg) for (e, a) in zip(self.engine_firing, self.acceleration_ms2)]
            self.location_km, self.velocity_ms = self.compute_simulation_step(self.location_km, self.velocity_ms, self.acceleration_ms2, timestep_s)
        if len(self.future) < 2:
            self.compute_future(timestep_s)
        self.location_km, self.velocity_ms = self.future[0]
        self.future = self.future[1:]
        last_location_km, last_velocity_ms = self.future[-1]
        last_acceleration_ms2 = self.planet.gravity_acceleration_components_ms2(last_location_km)
        new_location_km, new_velocity_ms = self.compute_simulation_step(last_location_km, last_velocity_ms, last_acceleration_ms2, timestep_s)
        self.future.append((new_location_km, new_velocity_ms))
        
    def compute_future(self, timestep_s = 10.0, depth = 15000):
        cur_location_km = self.location_km
        cur_velocity_ms = self.velocity_ms
        self.future = [(cur_location_km, cur_velocity_ms)]
        for _ in range(depth):
            cur_acceleration_ms2 = self.planet.gravity_acceleration_components_ms2(cur_location_km)
            cur_location_km, cur_velocity_ms = self.compute_simulation_step(cur_location_km, cur_velocity_ms, cur_acceleration_ms2, timestep_s)
            self.future.append((cur_location_km, cur_velocity_ms))

    def process_event(self, event):
        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_UP:
                self.engine_firing[1] = 1
            elif event.key == pygame.K_DOWN:
                self.engine_firing[1] = -1
            elif event.key == pygame.K_LEFT:
                self.engine_firing[0] = -1
            elif event.key == pygame.K_RIGHT:
                self.engine_firing[0] = 1
        elif event.type == pygame.KEYUP:
            if event.key in [pygame.K_UP, pygame.K_DOWN]:
                self.engine_firing[1] = 0
            elif event.key in [pygame.K_LEFT, pygame.K_RIGHT]:
                self.engine_firing[0] = 0


class Scene:
    def __init__(self, screen_size, fill_color=(255,255,255)):
        self.screen_size = screen_size
        self.layers = []
        self.fill_color = fill_color
        
    def add_layer(self, layer, align="TOPLEFT"):
        self.layers.append((layer, align))
    
    def set_screen_size(self, screen_size):
        self.screen_size = screen_size
        for (layer, align) in self.layers:
            layer.screen_size = self.screen_size
    
    def render(self, target_surface, pixels_per_km=1):
        for (layer, align) in self.layers:
            layer_size = layer.get_size()
            target_size = target_surface.get_size()
            
            #A hack to make scale work correctly
            large_target = pygame.Surface(layer_size)
            layer.render(large_target, (0,0), pixels_per_km)
            layer_size = layer.get_size()

            if align == "CENTER":
                dest = [int( (target_size[_x] - layer_size[_x])/2 ) for _x in [0, 1]]
            else:
                dest = [0, 0]
        
            if min(dest) < 0:
                #Crop the surface because pygame can't deal with negative dest
                large_target = pygame.Surface(layer_size)
                large_target.fill(self.fill_color)
                large_target.set_colorkey(self.fill_color)
                layer.render(large_target, (0,0), pixels_per_km)
                crop_dest = [-min(dest[_x], 0) for _x in [0, 1]]
                target_surface.blit(large_target, (max(dest[0],0), max(dest[1], 0)), area=(crop_dest[0], crop_dest[1], target_size[0], target_size[1]))
            else:
                layer.render(target_surface, dest, pixels_per_km)
    

class Simulation:
    def __init__(self, time_acceleration = 10):
        self.sputniks = []
        self.time_acceleration = time_acceleration
        
    def add_spacecraft(self, sputnik):
        self.sputniks.append(sputnik)
    
    def process_event(self, event):
        if len(self.sputniks) >= 1:
            self.sputniks[0].process_event(event)
    
    def make_simulation_step(self):
        for i in range(len(self.sputniks)):
            for _ in range(self.time_acceleration):
                self.sputniks[i].make_simulation_step(timestep_s = 1.0)
                if self.sputniks[i].engine_firing[0] != 0 or self.sputniks[i].engine_firing[1] != 0:
                    break
            


def main():
    #Start pygame
    pygame.init()
    
    
    #Get the current screen size, set default
    display_info = pygame.display.Info()
    default_scale = 0.75
    default_size = int(min(display_info.current_w, display_info.current_h)*default_scale)
    screen_size = [default_size, default_size]
    pixels_per_km=1.0/32.0


    #Create the main window
    main_surface = pygame.display.set_mode(screen_size, pygame.RESIZABLE)
    pygame.display.set_caption("Sputnik")
    colorkey = (255, 255, 255)
    main_surface.fill(colorkey)
    main_surface.set_colorkey(colorkey)
    
    #Create celestial objects
    planet = Planet()
    
    
    #Create spacecraft
    sputnik1 = Sputnik(orbiting_planet = planet, screen_size=screen_size, start_altitude_km = 3000, start_velocity_ms = 7000)
    
    
    #Create simulation
    simulation = Simulation()
    simulation.add_spacecraft(sputnik1)
    
    
    #Create the scene
    scene = Scene(screen_size, colorkey)
    background = Background(screen_size)
    scene.add_layer(background, align="CENTER")
    scene.add_layer(planet, align="CENTER")
    scene.add_layer(sputnik1, align="TOPLEFT")
    
    
    #Main loop
    while True:
        #Render the scene
        scene.render(main_surface, pixels_per_km)
        
        #Show everything to the user
        pygame.display.update()
    
        #Process major events
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            if event.type == pygame.VIDEORESIZE:
                screen_size = [event.w, event.h]
                scene.set_screen_size(screen_size)
                main_surface = pygame.display.set_mode(screen_size, pygame.RESIZABLE)
            if event.type == pygame.MOUSEWHEEL:
                pixels_per_km += pixels_per_km*event.y/10
            if event.type in [pygame.KEYDOWN, pygame.KEYUP] and event.key in [pygame.K_UP, pygame.K_DOWN, pygame.K_LEFT, pygame.K_RIGHT]:
                simulation.process_event(event)

        #Run the simulation
        simulation.make_simulation_step()
        
    #Just in case we somehow get to here, exit
    pygame.quit()
    sys.exit()

if __name__ == "__main__":
    main()