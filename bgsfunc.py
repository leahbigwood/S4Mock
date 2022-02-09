import numpy as np
from astropy.table import Table

def tile2rosette(tile):
    if tile < 433:
        return (tile-1)//27
    else:
        if tile >= 433 and tile < 436:
            return 13
        if tile >= 436 and tile < 439:
            return 14
        if tile >= 439 and tile < 442:
            return 15
        if tile >= 442 and tile <=480:
            return (tile-442)//3
            
        if tile > 480:
            return tile//30    
    return 999999 

def get_tiles(fpath='/global/cscratch1/sd/mjwilson/S4MOCK/tiles-sv3.ecsv'):
    
    tiles = Table.read(fpath)
    tiles = tiles[(tiles['STATUS'] == 'done') & (tiles['PROGRAM']=='BRIGHT')]
    tiles['ROSETTE'] = np.array([tile2rosette(x) for x in tiles['TILEID']])
    
    return tiles 

if __name__ == '__main__':
    # newfunc(2)
    
    tiles = get_tiles()
    
    tiles.print()
    