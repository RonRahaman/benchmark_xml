!!!! WRITE LUDEF, main
!!!! WRITE LUDEF, write_prolog
module xml_data_geometry_t
   use READ_XML_PRIMITIVES
   use WRITE_XML_PRIMITIVES
   use XMLPARSE
   implicit none
   save
   integer, private :: lurep_
   logical, private :: strict_
!!!! WRITE LUDEF, add_typedef

type cell_xml
!!!! WRITE LUDEF, add_variable
   integer                                         :: id
!!!! WRITE LUDEF, add_variable
   integer                                         :: universe
!!!! WRITE LUDEF, add_variable
   character(len=12)                                :: material
!!!! WRITE LUDEF, add_variable
   integer                                         :: fill
!!!! WRITE LUDEF, add_variable
   integer, dimension(:), pointer                  :: surfaces => null()
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0)), dimension(:), pointer   :: rotation => null()
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0)), dimension(:), pointer   :: translation => null()
end type cell_xml
!!!! WRITE LUDEF, add_typedef

type surface_xml
!!!! WRITE LUDEF, add_variable
   integer                                         :: id
!!!! WRITE LUDEF, add_variable
   character(len=15)                                :: type
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0)), dimension(:), pointer   :: coeffs => null()
!!!! WRITE LUDEF, add_variable
   character(len=12)                                :: boundary
end type surface_xml
!!!! WRITE LUDEF, add_typedef

type lattice_xml
!!!! WRITE LUDEF, add_variable
   integer                                         :: id
!!!! WRITE LUDEF, add_variable
   character(len=12)                                :: type
!!!! WRITE LUDEF, add_variable
   integer, dimension(:), pointer                  :: dimension => null()
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0)), dimension(:), pointer   :: lower_left => null()
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0)), dimension(:), pointer   :: width => null()
!!!! WRITE LUDEF, add_variable
   integer, dimension(:), pointer                  :: universes => null()
!!!! WRITE LUDEF, add_variable
   integer                                         :: outside
end type lattice_xml
!!!! WRITE LUDEF, add_variable
   type(cell_xml), dimension(:), pointer           :: cell_ => null()
!!!! WRITE LUDEF, add_variable
   type(surface_xml), dimension(:), pointer        :: surface_ => null()
!!!! WRITE LUDEF, add_variable
   type(lattice_xml), dimension(:), pointer        :: lattice_ => null()
!!!! WRITE LUSUBS, main
!!!! WRITE LUSUBS, write_prolog
contains
!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_cell_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar, count_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(cell_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar
   integer, intent(inout)                       :: count_dvar

   type(cell_xml), dimension(:), pointer :: temp_dvar

   count_dvar = count_dvar + 1
   do while (count_dvar .gt. size(dvar))
       allocate(temp_dvar(1:size(dvar)*2))
       temp_dvar(1:size(dvar)) = dvar
       deallocate(dvar)
       dvar => temp_dvar
       temp_dvar => null()
   enddo

   call read_xml_type_cell_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(count_dvar), has_dvar )
end subroutine read_xml_type_cell_xml_array

!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_cell_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(cell_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_id
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_universe
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_material
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_fill
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_surfaces
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_rotation
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_translation
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable
   has_id                               = .false.
!!!! WRITE LUSTART, add_variable
   has_universe                         = .false.
!!!! WRITE LUSTART, add_variable
   has_material                         = .false.
!!!! WRITE LUSTART, add_variable
   has_fill                             = .false.
!!!! WRITE LUSTART, add_variable
   has_surfaces                         = .false.
!!!! WRITE LUSTART, add_variable
   has_rotation                         = .false.
!!!! WRITE LUSTART, add_variable
   has_translation                      = .false.
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LULOOP, add_begin_loop
   call init_xml_type_cell_xml(dvar)
   has_dvar = .true.
!!!! WRITE LULOOP, add_begin_loop, component .eq. true
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
!!!! WRITE LULOOP, add_begin_loop, checktag .eq. true
      if ( endtag .and. tag == starttag ) then
         exit
      endif
!!!! WRITE LULOOP, add_begin_loop
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
!!!! WRITE LULOOP, add_variable
      case('id')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%id
!!!!        varname=id, vartype=integer, idx3= 5, idx6=-1, idx7=-1
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%id, has_id )
!!!! WRITE LULOOP, add_variable
      case('universe')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%universe
!!!!        varname=universe, vartype=integer, idx3= 5, idx6=-1, idx7=-1
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%universe, has_universe )
!!!! WRITE LULOOP, add_variable
      case('material')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%material
!!!!        varname=material, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%material, has_material )
!!!! WRITE LULOOP, add_variable
      case('fill')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%fill
!!!!        varname=fill, vartype=integer, idx3= 5, idx6=-1, idx7=-1
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%fill, has_fill )
!!!! WRITE LULOOP, add_variable
      case('surfaces')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%surfaces
!!!!        varname=surfaces, vartype=integer-array, idx3= 7, idx6=-1, idx7=-1
         call read_xml_integer_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%surfaces, has_surfaces )
!!!! WRITE LULOOP, add_variable
      case('rotation')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%rotation
!!!!        varname=rotation, vartype=double-array, idx3=15, idx6=-1, idx7=-1
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%rotation, has_rotation )
!!!! WRITE LULOOP, add_variable
      case('translation')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%translation
!!!!        varname=translation, vartype=double-array, idx3=15, idx6=-1, idx7=-1
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%translation, has_translation )
!!!! WRITE LUEND, open_tmp_files
!!!! WRITE LUEND, add_end_loop
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
!!!! WRITE LUEND, add_end_loop
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_id ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on id')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_surfaces ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on surfaces')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_rotation ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on rotation')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_translation ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on translation')
   endif
end subroutine read_xml_type_cell_xml
!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUDEFLT, add_typedef
subroutine init_xml_type_cell_xml_array( dvar )
   type(cell_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(1) )
end subroutine init_xml_type_cell_xml_array
subroutine init_xml_type_cell_xml(dvar)
   type(cell_xml) :: dvar
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%universe = 0
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%material = ''
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%fill = 0
end subroutine init_xml_type_cell_xml
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITE, add_typedef
subroutine write_xml_type_cell_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(cell_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_cell_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_cell_xml_array

subroutine write_xml_type_cell_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(cell_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_integer( info, 'id', indent+3, dvar%id)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_integer( info, 'universe', indent+3, dvar%universe)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'material', indent+3, dvar%material)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_integer( info, 'fill', indent+3, dvar%fill)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%surfaces)) call write_to_xml_integer_array( info, 'surfaces', indent+3, dvar%surfaces)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%rotation)) call write_to_xml_double_array( info, 'rotation', indent+3, dvar%rotation)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%translation)) call write_to_xml_double_array( info, 'translation', indent+3, dvar%translation)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_cell_xml

!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_surface_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar, count_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(surface_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar
   integer, intent(inout)                       :: count_dvar

   type(surface_xml), dimension(:), pointer :: temp_dvar

   count_dvar = count_dvar + 1
   do while (count_dvar .gt. size(dvar))
       allocate(temp_dvar(1:size(dvar)*2))
       temp_dvar(1:size(dvar)) = dvar
       deallocate(dvar)
       dvar => temp_dvar
       temp_dvar => null()
   enddo

   call read_xml_type_surface_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(count_dvar), has_dvar )
end subroutine read_xml_type_surface_xml_array

!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_surface_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(surface_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_id
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_type
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_coeffs
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_boundary
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable
   has_id                               = .false.
!!!! WRITE LUSTART, add_variable
   has_type                             = .false.
!!!! WRITE LUSTART, add_variable
   has_coeffs                           = .false.
!!!! WRITE LUSTART, add_variable
   has_boundary                         = .false.
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LULOOP, add_begin_loop
   call init_xml_type_surface_xml(dvar)
   has_dvar = .true.
!!!! WRITE LULOOP, add_begin_loop, component .eq. true
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
!!!! WRITE LULOOP, add_begin_loop, checktag .eq. true
      if ( endtag .and. tag == starttag ) then
         exit
      endif
!!!! WRITE LULOOP, add_begin_loop
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
!!!! WRITE LULOOP, add_variable
      case('id')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%id
!!!!        varname=id, vartype=integer, idx3= 5, idx6=-1, idx7=-1
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%id, has_id )
!!!! WRITE LULOOP, add_variable
      case('type')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%type
!!!!        varname=type, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type, has_type )
!!!! WRITE LULOOP, add_variable
      case('coeffs')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%coeffs
!!!!        varname=coeffs, vartype=double-array, idx3=15, idx6=-1, idx7=-1
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%coeffs, has_coeffs )
!!!! WRITE LULOOP, add_variable
      case('boundary')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%boundary
!!!!        varname=boundary, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%boundary, has_boundary )
!!!! WRITE LUEND, open_tmp_files
!!!! WRITE LUEND, add_end_loop
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
!!!! WRITE LUEND, add_end_loop
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_id ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on id')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_type ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on type')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_coeffs ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on coeffs')
   endif
end subroutine read_xml_type_surface_xml
!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUDEFLT, add_typedef
subroutine init_xml_type_surface_xml_array( dvar )
   type(surface_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(1) )
end subroutine init_xml_type_surface_xml_array
subroutine init_xml_type_surface_xml(dvar)
   type(surface_xml) :: dvar
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%boundary = 'transmit'
end subroutine init_xml_type_surface_xml
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITE, add_typedef
subroutine write_xml_type_surface_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(surface_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_surface_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_surface_xml_array

subroutine write_xml_type_surface_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(surface_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_integer( info, 'id', indent+3, dvar%id)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'type', indent+3, dvar%type)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%coeffs)) call write_to_xml_double_array( info, 'coeffs', indent+3, dvar%coeffs)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'boundary', indent+3, dvar%boundary)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_surface_xml

!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_lattice_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar, count_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(lattice_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar
   integer, intent(inout)                       :: count_dvar

   type(lattice_xml), dimension(:), pointer :: temp_dvar

   count_dvar = count_dvar + 1
   do while (count_dvar .gt. size(dvar))
       allocate(temp_dvar(1:size(dvar)*2))
       temp_dvar(1:size(dvar)) = dvar
       deallocate(dvar)
       dvar => temp_dvar
       temp_dvar => null()
   enddo

   call read_xml_type_lattice_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(count_dvar), has_dvar )
end subroutine read_xml_type_lattice_xml_array

!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_lattice_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(lattice_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_id
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_type
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_dimension
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_lower_left
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_width
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_universes
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_outside
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable
   has_id                               = .false.
!!!! WRITE LUSTART, add_variable
   has_type                             = .false.
!!!! WRITE LUSTART, add_variable
   has_dimension                        = .false.
!!!! WRITE LUSTART, add_variable
   has_lower_left                       = .false.
!!!! WRITE LUSTART, add_variable
   has_width                            = .false.
!!!! WRITE LUSTART, add_variable
   has_universes                        = .false.
!!!! WRITE LUSTART, add_variable
   has_outside                          = .false.
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LULOOP, add_begin_loop
   call init_xml_type_lattice_xml(dvar)
   has_dvar = .true.
!!!! WRITE LULOOP, add_begin_loop, component .eq. true
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
!!!! WRITE LULOOP, add_begin_loop, checktag .eq. true
      if ( endtag .and. tag == starttag ) then
         exit
      endif
!!!! WRITE LULOOP, add_begin_loop
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
!!!! WRITE LULOOP, add_variable
      case('id')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%id
!!!!        varname=id, vartype=integer, idx3= 5, idx6=-1, idx7=-1
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%id, has_id )
!!!! WRITE LULOOP, add_variable
      case('type')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%type
!!!!        varname=type, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type, has_type )
!!!! WRITE LULOOP, add_variable
      case('dimension')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%dimension
!!!!        varname=dimension, vartype=integer-array, idx3= 7, idx6=-1, idx7=-1
         call read_xml_integer_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%dimension, has_dimension )
!!!! WRITE LULOOP, add_variable
      case('lower_left')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%lower_left
!!!!        varname=lower_left, vartype=double-array, idx3=15, idx6=-1, idx7=-1
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%lower_left, has_lower_left )
!!!! WRITE LULOOP, add_variable
      case('width')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%width
!!!!        varname=width, vartype=double-array, idx3=15, idx6=-1, idx7=-1
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%width, has_width )
!!!! WRITE LULOOP, add_variable
      case('universes')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%universes
!!!!        varname=universes, vartype=integer-array, idx3= 7, idx6=-1, idx7=-1
         call read_xml_integer_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%universes, has_universes )
!!!! WRITE LULOOP, add_variable
      case('outside')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%outside
!!!!        varname=outside, vartype=integer, idx3= 5, idx6=-1, idx7=-1
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%outside, has_outside )
!!!! WRITE LUEND, open_tmp_files
!!!! WRITE LUEND, add_end_loop
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
!!!! WRITE LUEND, add_end_loop
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_id ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on id')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_dimension ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on dimension')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_lower_left ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on lower_left')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_width ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on width')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_universes ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on universes')
   endif
end subroutine read_xml_type_lattice_xml
!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUDEFLT, add_typedef
subroutine init_xml_type_lattice_xml_array( dvar )
   type(lattice_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(1) )
end subroutine init_xml_type_lattice_xml_array
subroutine init_xml_type_lattice_xml(dvar)
   type(lattice_xml) :: dvar
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%type = 'rectangular'
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%outside = 0
end subroutine init_xml_type_lattice_xml
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITE, add_typedef
subroutine write_xml_type_lattice_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(lattice_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_lattice_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_lattice_xml_array

subroutine write_xml_type_lattice_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(lattice_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_integer( info, 'id', indent+3, dvar%id)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'type', indent+3, dvar%type)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%dimension)) call write_to_xml_integer_array( info, 'dimension', indent+3, dvar%dimension)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%lower_left)) call write_to_xml_double_array( info, 'lower_left', indent+3, dvar%lower_left)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%width)) call write_to_xml_double_array( info, 'width', indent+3, dvar%width)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%universes)) call write_to_xml_integer_array( info, 'universes', indent+3, dvar%universes)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_integer( info, 'outside', indent+3, dvar%outside)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_lattice_xml

!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, write_prolog
subroutine read_xml_file_geometry_t(fname, lurep, errout)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep
   logical, intent(out), optional         :: errout

   type(XML_PARSE)                        :: info
   logical                                :: error
   character(len=80)                      :: tag
   character(len=80)                      :: starttag
   logical                                :: endtag
   character(len=250), dimension(1:2,1:20) :: attribs
   integer                                :: noattribs
   character(len=1000), dimension(1:100000)  :: data
   integer                                :: nodata
!!!! WRITE LUPROLOG, add_variable, idx7>=1
   type(cell_xml), dimension(:), pointer :: temp_cell_ => null()
   integer :: count_cell_
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_cell_
!!!! WRITE LUPROLOG, add_variable, idx7>=1
   type(surface_xml), dimension(:), pointer :: temp_surface_ => null()
   integer :: count_surface_
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_surface_
!!!! WRITE LUPROLOG, add_variable, idx7>=1
   type(lattice_xml), dimension(:), pointer :: temp_lattice_ => null()
   integer :: count_lattice_
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_lattice_
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_cell_ = 0
!!!! WRITE LUSTART, add_variable
   has_cell_                            = .false.
!!!! WRITE LUSTART, add_variable
   allocate(cell_(1))
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_surface_ = 0
!!!! WRITE LUSTART, add_variable
   has_surface_                         = .false.
!!!! WRITE LUSTART, add_variable
   allocate(surface_(1))
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_lattice_ = 0
!!!! WRITE LUSTART, add_variable
   has_lattice_                         = .false.
!!!! WRITE LUSTART, add_variable
   allocate(lattice_(1))
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LUMAIN, main
!!!! WRITE LUMAIN, write_prolog

   call init_xml_file_geometry_t
   call xml_open( info, fname, .true. )
   call xml_options( info, report_errors=.false., ignore_whitespace=.true.)
   lurep_ = 0
   if ( present(lurep) ) then
      lurep_ = lurep
      call xml_options( info, report_lun=lurep )
   endif
   do
      call xml_get( info, starttag, endtag, attribs, noattribs, &
         data, nodata)
      if ( starttag /= '!--' ) exit
   enddo
   if ( starttag /= "geometry" ) then
      call xml_report_errors( info, &
         'XML-file should have root element "geometry"')
      error = .true.
      call xml_close(info)
      return
   endif
   strict_ = .false.
!!!! WRITE LULOOP, add_begin_loop, component .eq. false
   error = .false.
   do
      call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
      if ( xml_error(info) ) then
         write(lurep_,*) 'Error reading input file!'
         error = .true.
         return
      endif
!!!! WRITE LULOOP, add_begin_loop, checktag .eq. true
      if ( endtag .and. tag == starttag ) then
         exit
      endif
!!!! WRITE LULOOP, add_begin_loop
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
!!!! WRITE LULOOP, add_variable
      case('cell')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=cell_
!!!!        varname=cell_, vartype=cell_xml-array, idx3=29, idx6=-1, idx7= 4
         call read_xml_type_cell_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            cell_, has_cell_, count_cell_ )
!!!! WRITE LULOOP, add_variable
      case('surface')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=surface_
!!!!        varname=surface_, vartype=surface_xml-array, idx3=32, idx6=-1, idx7= 4
         call read_xml_type_surface_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            surface_, has_surface_, count_surface_ )
!!!! WRITE LULOOP, add_variable
      case('lattice')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=lattice_
!!!!        varname=lattice_, vartype=lattice_xml-array, idx3=35, idx6=-1, idx7= 4
         call read_xml_type_lattice_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            lattice_, has_lattice_, count_lattice_ )
!!!! WRITE LUEND, open_tmp_files
!!!! WRITE LUEND, add_end_loop
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
!!!! WRITE LUEND, add_end_loop
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_cell_ ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=false
      error = .true.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on cell_')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_cell_(1:count_cell_))
   temp_cell_ = cell_(1:count_cell_)
   deallocate(cell_)
   cell_ => temp_cell_
   temp_cell_ => null()
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_surface_ ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=false
      error = .true.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on surface_')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_surface_(1:count_surface_))
   temp_surface_ = surface_(1:count_surface_)
   deallocate(surface_)
   surface_ => temp_surface_
   temp_surface_ => null()
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_lattice_ ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=false
      error = .true.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on lattice_')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_lattice_(1:count_lattice_))
   temp_lattice_ = lattice_(1:count_lattice_)
   deallocate(lattice_)
   lattice_ => temp_lattice_
   temp_lattice_ => null()
   if ( present(errout) ) errout = error
   call xml_close(info)
end subroutine

!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITP, main
!!!! WRITE LUWRV, main
!!!! WRITE LUWRITP, write_prolog
subroutine write_xml_file_geometry_t(fname, lurep)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep

   type(XML_PARSE)                        :: info
   integer                                :: indent = 0

   call xml_open( info, fname, .false. )
   call xml_options( info, report_errors=.true.)
   if ( present(lurep) ) then
       call xml_options( info, report_errors=.true.)
   endif
   write(info%lun,'(a)') &
      '<geometry>'
!!!! WRITE LUWRV, add_variable, component=false
   if(associated(cell_))   call write_xml_type_cell_xml_array( info, 'cell', indent+3, cell_)
!!!! WRITE LUWRV, add_variable, component=false
   if(associated(surface_))   call write_xml_type_surface_xml_array( info, 'surface', indent+3, surface_)
!!!! WRITE LUWRV, add_variable, component=false
   if(associated(lattice_))   call write_xml_type_lattice_xml_array( info, 'lattice', indent+3, lattice_)
   write(info%lun,'(a)') '</geometry>'
   call xml_close(info)
end subroutine

!!!! WRITE LUINIT, main
subroutine init_xml_file_geometry_t

end subroutine

end module
