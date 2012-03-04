# read in Region 1 description
gfx read node 2DCoupled-LinearElasticity-LinearElasticity_1.part0.exnode
gfx read element 2DCoupled-LinearElasticity-LinearElasticity_1.part0.exelem

# define deformed geometry
gfx define field "deformed_geom1" component Dependent1.1 Dependent1.2

# display deformed geometry
gfx define faces egroup "Region1"
gfx modify g_element "Region1" lines coordinate deformed_geom1 select_on material default selected_material default_selected
#gfx modify g_element "Region1" surfaces coordinate deformed_geom1 select_on material tissue selected_material default_selected render_shaded
gfx modify g_element "Region1" node_points coordinate deformed_geom1 glyph sphere General size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default selected_material default_selected

# display undeformed lines
gfx modify g_element "Region1" lines select_on material green selected_material default_selected

# read in Region 2 description
gfx read node 2DCoupled-LinearElasticity-LinearElasticity_2.part0.exnode node_offset 1000
gfx read element 2DCoupled-LinearElasticity-LinearElasticity_2.part0.exelem node_offset 1000 element_offset 1000 line_offset 1000 face_offset 1000

# define deformed geometry
gfx define field "deformed_geom2" component Dependent2.1 Dependent2.2

# display deformed geometry
gfx define faces egroup "Region2"
gfx modify g_element "Region2" lines coordinate deformed_geom2 select_on material default selected_material default_selected
#gfx modify g_element "Region2" surfaces coordinate deformed_geom2 select_on material tissue selected_material default_selected render_shaded
gfx modify g_element "Region2" node_points coordinate deformed_geom2 glyph sphere General size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default selected_material default_selected

# display undeformed lines
gfx modify g_element "Region2" lines select_on material green selected_material default_selected

gfx create window 1

gfx create axes length 5 material default
gfx draw axes

gfx edit scene
gfx modify window 1 set antialias 2
