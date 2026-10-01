module mesh_module
  use precision_module
  implicit none

  private

  public :: mesh_type

  public :: vert_type
  public :: face_type
  public :: minimal_face_type
  public :: boundary_face_type
  public :: sub_face_type
  public :: elem_type
  public :: sub_elem_type

  type :: vert_type
    logical :: is_bound = .FALSE., is_ghost = .FALSE.
    integer(kind=ENTIER) :: n_sub_faces_neigh, n_sub_elems_neigh
    integer(kind=ENTIER) :: n_elems_neigh, n_faces_neigh
    integer(kind=ENTIER) :: id_glob
    integer(kind=ENTIER), dimension(:), allocatable :: sub_face_neigh, sub_elem_neigh
    integer(kind=ENTIER), dimension(:), allocatable :: elem_neigh, face_neigh
    real(kind=DOUBLE), dimension(3) :: coord
    real(kind=DOUBLE) :: volume
  end type vert_type

  type :: face_type
    integer(kind=ENTIER) :: n_vert = -1, n_sub_faces = -1
    integer(kind=ENTIER) :: left_neigh = 0, left_neigh_face = 0
    integer(kind=ENTIER) :: right_neigh = 0, right_neigh_face = 0
    integer(kind=ENTIER), dimension(:), allocatable :: sub_face, vert
    real(kind=DOUBLE) :: area
    real(kind=DOUBLE), dimension(3) :: coord, norm
  end type face_type

  type :: minimal_face_type
    integer(kind=ENTIER) :: n_vert = -1
    integer(kind=ENTIER) :: left_neigh = 0, left_neigh_face = 0
    integer(kind=ENTIER) :: right_neigh = 0, right_neigh_face = 0
    integer(kind=ENTIER), dimension(:), allocatable :: vert
  end type minimal_face_type

  type :: boundary_face_type
    integer(kind=ENTIER) :: n_vert = -1
    integer(kind=ENTIER), dimension(:), allocatable :: vert
    integer(kind=ENTIER) :: tag
  end type boundary_face_type

  type :: sub_face_type
    integer(kind=ENTIER) :: mesh_vert, mesh_face
    integer(kind=ENTIER) :: id_loc_around_node
    integer(kind=ENTIER) :: left_elem_neigh = 0, right_elem_neigh = 0
    integer(kind=ENTIER) :: left_sub_elem_neigh = 0, right_sub_elem_neigh = 0
    real(kind=DOUBLE) :: area
    real(kind=DOUBLE), dimension(3) :: norm
  end type sub_face_type

  type :: elem_type
    integer(kind=ENTIER) :: n_vert, n_faces, n_sub_faces, n_sub_elems, elem_kind, tag = -1
    integer(kind=ENTIER) :: n_neigh_by_vert
    integer(kind=ENTIER), dimension(:), allocatable :: vert, face, sub_face, sub_elem
    integer(kind=ENTIER), dimension(:), allocatable :: neigh_by_vert
    real(kind=DOUBLE) :: volume
    real(kind=DOUBLE), dimension(3) :: coord
    logical :: is_ghost = .FALSE.
  end type elem_type

  type :: sub_elem_type
    integer(kind=ENTIER) :: id_loc_around_node
    integer(kind=ENTIER) :: mesh_vert, mesh_elem, n_sub_faces
    integer(kind=ENTIER), dimension(:), allocatable :: sub_face
    real(kind=DOUBLE) :: volume
  end type sub_elem_type

  type :: mesh_type
    integer(kind=ENTIER) :: n_vert, n_elems, n_faces, n_interior_elems = 0
    integer(kind=ENTIER) :: n_sub_faces, n_boundary_faces, n_sub_elems
    integer(kind=ENTIER) :: n_max_sub_faces_around_node
    type(vert_type), dimension(:), allocatable :: vert
    type(elem_type), dimension(:), allocatable :: elem
    type(face_type), dimension(:), allocatable :: face
    type(sub_face_type), dimension(:), allocatable :: sub_face
    type(boundary_face_type), dimension(:), allocatable :: boundary_face
    type(sub_elem_type), dimension(:), allocatable :: sub_elem

    !For reading only, gets transfered to face
    type(minimal_face_type), dimension(:), allocatable :: minimal_face
  end type mesh_type
end module mesh_module
