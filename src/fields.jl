"""
    Cell

A type describing the location at the center of a grid cell.
"""
struct Cell end

"""
	Face

A type describing the location at the face of a grid cell.
"""
struct Face end

"""
    Field{LX, LY, LZ, A, G} <: AbstractField{A, G}

A field defined at the location (`LX`, `LY`, `LZ`) which can be either `Cell` or `Face`.
"""
struct Field{Lx, Ly, Lz, A, G} <: AbstractField{A, G}
    data :: A
    grid :: G
end

"""
	CellField

A field defined at the cell centers. Used for pressure and tracers.
"""
const CellField  = Field{Cell, Cell, Cell}

"""
	FaceFieldX

A field defined at the faces along the x-direction. Used for horizontal velocity u.
"""
const FaceFieldX = Field{Face, Cell, Cell}

"""
	FaceFieldY

A field defined at the faces along the y-direction. Used for horizontal velocity v.
"""
const FaceFieldY = Field{Cell, Face, Cell}

"""
	FaceFieldZ

A field defined at the faces along the z-direction. Used for vertical velocity w.
"""
const FaceFieldZ = Field{Cell, Cell, Face}

