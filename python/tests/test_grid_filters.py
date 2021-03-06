import pytest
import numpy as np

from damask import grid_filters

class TestGridFilters:

    def test_cell_coord0(self):
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         coord = grid_filters.cell_coord0(grid,size)
         assert np.allclose(coord[0,0,0],size/grid*.5) and coord.shape == tuple(grid) + (3,)

    def test_node_coord0(self):
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         coord = grid_filters.node_coord0(grid,size)
         assert np.allclose(coord[-1,-1,-1],size) and coord.shape == tuple(grid+1) + (3,)

    def test_coord0(self):
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         c = grid_filters.cell_coord0(grid+1,size+size/grid)
         n = grid_filters.node_coord0(grid,size) + size/grid*.5
         assert np.allclose(c,n)

    @pytest.mark.parametrize('mode',['cell','node'])
    def test_grid_DNA(self,mode):
         """Ensure that xx_coord0_gridSizeOrigin is the inverse of xx_coord0."""
         grid   = np.random.randint(8,32,(3))
         size   = np.random.random(3)
         origin = np.random.random(3)
         coord0 = eval(f'grid_filters.{mode}_coord0(grid,size,origin)')                     # noqa
         _grid,_size,_origin = eval(f'grid_filters.{mode}_coord0_gridSizeOrigin(coord0.reshape(-1,3,order="F"))')
         assert np.allclose(grid,_grid) and np.allclose(size,_size) and np.allclose(origin,_origin)

    def test_displacement_fluct_equivalence(self):
         """Ensure that fluctuations are periodic."""
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         F    = np.random.random(tuple(grid)+(3,3))
         assert np.allclose(grid_filters.node_displacement_fluct(size,F),
                            grid_filters.cell_2_node(grid_filters.cell_displacement_fluct(size,F)))

    def test_interpolation_to_node(self):
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         F    = np.random.random(tuple(grid)+(3,3))
         assert np.allclose(grid_filters.node_coord(size,F) [1:-1,1:-1,1:-1],
                            grid_filters.cell_2_node(grid_filters.cell_coord(size,F))[1:-1,1:-1,1:-1])

    def test_interpolation_to_cell(self):
         grid = np.random.randint(1,30,(3))

         node_coord_x = np.linspace(0,np.pi*2,num=grid[0]+1)
         node_field_x = np.cos(node_coord_x)
         node_field   = np.broadcast_to(node_field_x.reshape(-1,1,1),grid+1)

         cell_coord_x = node_coord_x[:-1]+node_coord_x[1]*.5
         cell_field_x = np.interp(cell_coord_x,node_coord_x,node_field_x,period=np.pi*2.)
         cell_field   = np.broadcast_to(cell_field_x.reshape(-1,1,1),grid)

         assert np.allclose(cell_field,grid_filters.node_2_cell(node_field))

    @pytest.mark.parametrize('mode',['cell','node'])
    def test_coord0_origin(self,mode):
         origin= np.random.random(3)
         size  = np.random.random(3)                                                                # noqa
         grid  = np.random.randint(8,32,(3))
         shifted   = eval(f'grid_filters.{mode}_coord0(grid,size,origin)')
         unshifted = eval(f'grid_filters.{mode}_coord0(grid,size)')
         if   mode == 'cell':
            assert  np.allclose(shifted,unshifted+np.broadcast_to(origin,tuple(grid)  +(3,)))
         elif mode == 'node':
            assert  np.allclose(shifted,unshifted+np.broadcast_to(origin,tuple(grid+1)+(3,)))

    @pytest.mark.parametrize('function',[grid_filters.cell_displacement_avg,
                                         grid_filters.node_displacement_avg])
    def test_displacement_avg_vanishes(self,function):
         """Ensure that random fluctuations in F do not result in average displacement."""
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         F    = np.random.random(tuple(grid)+(3,3))
         F   += np.eye(3) - np.average(F,axis=(0,1,2))
         assert np.allclose(function(size,F),0.0)

    @pytest.mark.parametrize('function',[grid_filters.cell_displacement_fluct,
                                         grid_filters.node_displacement_fluct])
    def test_displacement_fluct_vanishes(self,function):
         """Ensure that constant F does not result in fluctuating displacement."""
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         F    = np.broadcast_to(np.random.random((3,3)), tuple(grid)+(3,3))
         assert np.allclose(function(size,F),0.0)

    @pytest.mark.parametrize('function',[grid_filters.coord0_check,
                                         grid_filters.node_coord0_gridSizeOrigin,
                                         grid_filters.cell_coord0_gridSizeOrigin])
    def test_invalid_coordinates(self,function):
        invalid_coordinates = np.random.random((np.random.randint(12,52),3))
        with pytest.raises(ValueError):
            function(invalid_coordinates)

    @pytest.mark.parametrize('function',[grid_filters.node_coord0_gridSizeOrigin,
                                         grid_filters.cell_coord0_gridSizeOrigin])
    def test_uneven_spaced_coordinates(self,function):
        start = np.random.random(3)
        end   = np.random.random(3)*10. + start
        grid  = np.random.randint(8,32,(3))
        uneven = np.stack(np.meshgrid(np.logspace(start[0],end[0],grid[0]),
                                      np.logspace(start[1],end[1],grid[1]),
                                      np.logspace(start[2],end[2],grid[2]),indexing = 'ij'),
                           axis = -1).reshape((grid.prod(),3),order='F')
        with pytest.raises(ValueError):
            function(uneven)


    @pytest.mark.parametrize('mode',[True,False])
    @pytest.mark.parametrize('function',[grid_filters.node_coord0_gridSizeOrigin,
                                         grid_filters.cell_coord0_gridSizeOrigin])
    def test_unordered_coordinates(self,function,mode):
        origin = np.random.random(3)
        size   = np.random.random(3)*10.+origin
        grid  = np.random.randint(8,32,(3))
        unordered = grid_filters.node_coord0(grid,size,origin).reshape(-1,3)
        if mode:
            with pytest.raises(ValueError):
                function(unordered,mode)
        else:
            function(unordered,mode)

    def test_regrid(self):
         size = np.random.random(3)
         grid = np.random.randint(8,32,(3))
         F    = np.broadcast_to(np.eye(3), tuple(grid)+(3,3))
         assert all(grid_filters.regrid(size,F,grid) == np.arange(grid.prod()))


    @pytest.mark.parametrize('differential_operator',[grid_filters.curl,
                                                      grid_filters.divergence,
                                                      grid_filters.gradient])
    def test_differential_operator_constant(self,differential_operator):
        size = np.random.random(3)+1.0
        grid = np.random.randint(8,32,(3))
        shapes = {
                  grid_filters.curl:      [(3,),(3,3)],
                  grid_filters.divergence:[(3,),(3,3)],
                  grid_filters.gradient:  [(1,),(3,)]
                 }
        for shape in shapes[differential_operator]:
            field = np.ones(tuple(grid)+shape)*np.random.random()*1.0e5
            assert np.allclose(differential_operator(size,field),0.0)


    grad_test_data = [
    (['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0'],
     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0', '0.0',
      '0.0',                                                  '0.0', '0.0',
      '0.0',                                                  '0.0', '0.0']),

    (['0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0' ],
     ['0.0', '0.0',                                                   '0.0',
      '0.0', '-np.pi*2/size[1]*np.sin(np.pi*2*nodes[...,1]/size[1])', '0.0',
      '0.0', '0.0',                                                   '0.0' ]),

    (['1.0', '0.0', '2.0*np.cos(np.pi*2*nodes[...,2]/size[2])'],
     ['0.0', '0.0', '0.0',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', '-2.0*np.pi*2/size[2]*np.sin(np.pi*2*nodes[...,2]/size[2])']),

    (['np.cos(np.pi*2*nodes[...,2]/size[2])', '3.0', 'np.sin(np.pi*2*nodes[...,2]/size[2])'],
     ['0.0', '0.0', '-np.sin(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', ' np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']),

    (['np.sin(np.pi*2*nodes[...,0]/size[0])',
      'np.sin(np.pi*2*nodes[...,1]/size[1])',
      'np.sin(np.pi*2*nodes[...,2]/size[2])'],
     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0', '0.0',
      '0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]', '0.0',
      '0.0', '0.0', 'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']),

    (['np.sin(np.pi*2*nodes[...,0]/size[0])'],
     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]', '0.0', '0.0']),

    (['8.0'],
     ['0.0', '0.0', '0.0' ])
                    ]

    @pytest.mark.parametrize('field_def,grad_def',grad_test_data)
    def test_grad(self,field_def,grad_def):
        size = np.random.random(3)+1.0
        grid = np.random.randint(8,32,(3))

        nodes = grid_filters.cell_coord0(grid,size)
        my_locals = locals()                                                                        # needed for list comprehension

        field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),grid) for f in field_def],axis=-1)
        field = field.reshape(tuple(grid) + ((3,) if len(field_def)==3 else (1,)))
        grad = np.stack([np.broadcast_to(eval(c,globals(),my_locals),grid) for c in grad_def], axis=-1)
        grad = grad.reshape(tuple(grid) + ((3,3) if len(grad_def)==9 else (3,)))

        assert np.allclose(grad,grid_filters.gradient(size,field))


    curl_test_data = [
    (['np.sin(np.pi*2*nodes[...,2]/size[2])', '0.0', '0.0',
      '0.0',                                  '0.0', '0.0',
      '0.0',                                  '0.0', '0.0'],
     ['0.0'                                                 , '0.0', '0.0',
      'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]', '0.0', '0.0',
      '0.0',                                                  '0.0', '0.0']),

    (['np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0', '0.0',
      '0.0',                                  '0.0', '0.0',
      'np.cos(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0'],
     ['0.0',                                                  '0.0', '0.0',
      '0.0',                                                  '0.0', '0.0',
      'np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]', '0.0', '0.0']),

    (['np.sin(np.pi*2*nodes[...,0]/size[0])','np.cos(np.pi*2*nodes[...,1]/size[1])','np.sin(np.pi*2*nodes[...,2]/size[2])',
      'np.sin(np.pi*2*nodes[...,0]/size[0])','np.cos(np.pi*2*nodes[...,1]/size[1])','np.sin(np.pi*2*nodes[...,2]/size[2])',
      'np.sin(np.pi*2*nodes[...,0]/size[0])','np.cos(np.pi*2*nodes[...,1]/size[1])','np.sin(np.pi*2*nodes[...,2]/size[2])'],
     ['0.0', '0.0', '0.0',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', '0.0']),

    (['5.0', '0.0', '0.0',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', '2*np.cos(np.pi*2*nodes[...,1]/size[1])'],
     ['0.0', '0.0', '-2*np.pi*2/size[1]*np.sin(np.pi*2*nodes[...,1]/size[1])',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', '0.0']),

    ([ '4*np.sin(np.pi*2*nodes[...,2]/size[2])',
       '8*np.sin(np.pi*2*nodes[...,0]/size[0])',
      '16*np.sin(np.pi*2*nodes[...,1]/size[1])'],
     ['16*np.pi*2/size[1]*np.cos(np.pi*2*nodes[...,1]/size[1])',
       '4*np.pi*2/size[2]*np.cos(np.pi*2*nodes[...,2]/size[2])',
       '8*np.pi*2/size[0]*np.cos(np.pi*2*nodes[...,0]/size[0])']),

    (['0.0',
      'np.cos(np.pi*2*nodes[...,0]/size[0])+5*np.cos(np.pi*2*nodes[...,2]/size[2])',
      '0.0'],
     ['5*np.sin(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',
      '0.0',
      '-np.sin(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]'])
                     ]

    @pytest.mark.parametrize('field_def,curl_def',curl_test_data)
    def test_curl(self,field_def,curl_def):
        size = np.random.random(3)+1.0
        grid = np.random.randint(8,32,(3))

        nodes = grid_filters.cell_coord0(grid,size)
        my_locals = locals()                                                                        # needed for list comprehension

        field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),grid) for f in field_def],axis=-1)
        field = field.reshape(tuple(grid) + ((3,3) if len(field_def)==9 else (3,)))
        curl = np.stack([np.broadcast_to(eval(c,globals(),my_locals),grid) for c in curl_def], axis=-1)
        curl = curl.reshape(tuple(grid) + ((3,3) if len(curl_def)==9 else (3,)))

        assert np.allclose(curl,grid_filters.curl(size,field))


    div_test_data =[
    (['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0',
      '0.0'                                 , '0.0', '0.0',
      '0.0'                                 , '0.0', '0.0'],
     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]','0.0', '0.0']),

    (['0.0', '0.0',                                  '0.0',
      '0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0',
      '0.0', '0.0',                                  '0.0'],
     ['0.0', '-np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]', '0.0']),

    (['1.0', '0.0', '0.0',
      '0.0', '0.0', '0.0',
      '0.0', '0.0', '2*np.cos(np.pi*2*nodes[...,2]/size[2])' ],
     ['0.0', '0.0', '-2.0*np.pi*2/size[2]*np.sin(np.pi*2*nodes[...,2]/size[2])']
     ),

    ([ '23.0', '0.0',   'np.sin(np.pi*2*nodes[...,2]/size[2])',
       '0.0',  '100.0', 'np.sin(np.pi*2*nodes[...,2]/size[2])',
       '0.0',  '0.0',   'np.sin(np.pi*2*nodes[...,2]/size[2])'],
      ['np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]',\
       'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]', \
       'np.cos(np.pi*2*nodes[...,2]/size[2])*np.pi*2/size[2]']),

    (['400.0',                                '0.0',                                  '0.0',
      'np.sin(np.pi*2*nodes[...,0]/size[0])', 'np.sin(np.pi*2*nodes[...,1]/size[1])', 'np.sin(np.pi*2*nodes[...,2]/size[2])',
      '0.0',                                  '10.0',                                 '6.0'],
     ['0.0','np.sum(np.cos(np.pi*2*nodes/size)*np.pi*2/size,axis=-1)', '0.0' ]),

    (['np.sin(np.pi*2*nodes[...,0]/size[0])', '0.0', '0.0'],
     ['np.cos(np.pi*2*nodes[...,0]/size[0])*np.pi*2/size[0]',]),

    (['0.0', 'np.cos(np.pi*2*nodes[...,1]/size[1])', '0.0' ],
     ['-np.sin(np.pi*2*nodes[...,1]/size[1])*np.pi*2/size[1]'])
     ]

    @pytest.mark.parametrize('field_def,div_def',div_test_data)

    def test_div(self,field_def,div_def):
        size = np.random.random(3)+1.0
        grid = np.random.randint(8,32,(3))

        nodes = grid_filters.cell_coord0(grid,size)
        my_locals = locals()                                                                        # needed for list comprehension

        field = np.stack([np.broadcast_to(eval(f,globals(),my_locals),grid) for f in field_def],axis=-1)
        field = field.reshape(tuple(grid) + ((3,3) if len(field_def)==9 else (3,)))
        div = np.stack([np.broadcast_to(eval(c,globals(),my_locals),grid) for c in div_def], axis=-1)
        if len(div_def)==3:
            div = div.reshape(tuple(grid) + ((3,)))
        else:
            div=div.reshape(tuple(grid))

        assert np.allclose(div,grid_filters.divergence(size,field))
