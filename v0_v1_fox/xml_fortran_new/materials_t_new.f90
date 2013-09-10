!!!! WRITE LUDEF, main
!!!! WRITE LUDEF, write_prolog
module xml_data_materials_t_new
   use READ_XML_PRIMITIVES
   use WRITE_XML_PRIMITIVES
   use XMLPARSE
   implicit none
   save
   integer, private :: lurep_
   logical, private :: strict_
!!!! WRITE LUDEF, add_typedef

type density_xml
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0))                          :: value
!!!! WRITE LUDEF, add_variable
   character(len=10)                                :: units
end type density_xml
!!!! WRITE LUDEF, add_typedef

type nuclide_xml
!!!! WRITE LUDEF, add_variable
   character(len=7)                                :: name
!!!! WRITE LUDEF, add_variable
   character(len=3)                                :: xs
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0))                          :: ao
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0))                          :: wo
end type nuclide_xml
!!!! WRITE LUDEF, add_typedef

type sab_xml
!!!! WRITE LUDEF, add_variable
   character(len=10)                                :: name
!!!! WRITE LUDEF, add_variable
   character(len=3)                                :: xs
end type sab_xml
!!!! WRITE LUDEF, add_typedef

type material_xml
!!!! WRITE LUDEF, add_variable
   integer                                         :: id
!!!! WRITE LUDEF, add_variable
   type(density_xml)                               :: density
!!!! WRITE LUDEF, add_variable
   type(nuclide_xml), dimension(:), pointer        :: nuclides => null()
!!!! WRITE LUDEF, add_variable
   type(nuclide_xml), dimension(:), pointer        :: elements => null()
!!!! WRITE LUDEF, add_variable
   type(sab_xml), dimension(:), pointer            :: sab => null()
end type material_xml
!!!! WRITE LUDEF, add_variable
   type(material_xml), dimension(:), pointer       :: material_ => null()
!!!! WRITE LUDEF, add_variable
   character(len=3)                                :: default_xs_
!!!! WRITE LUSUBS, main
!!!! WRITE LUSUBS, write_prolog
contains
!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_density_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar, count_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(density_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar
   integer, intent(inout)                       :: count_dvar

   integer                                      :: newsize
   integer                                      :: size_dvar
   type(density_xml), dimension(:), pointer :: newvar

   count_dvar = count_dvar + 1
   size_dvar = size(dvar)
   if (count_dvar .gt. size_dvar) then
       newsize = size_dvar * 2
       allocate(newvar(1:newsize))
       newvar(1:size_dvar) = dvar
       deallocate(dvar)
       dvar => newvar
   endif

   call read_xml_type_density_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(count_dvar), has_dvar )
end subroutine read_xml_type_density_xml_array

!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_density_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(density_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_value
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_units
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable
   has_value                            = .false.
!!!! WRITE LUSTART, add_variable
   has_units                            = .false.
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LULOOP, add_begin_loop
   call init_xml_type_density_xml(dvar)
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
      case('value')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%value
!!!!        varname=value, vartype=double, idx3=13, idx6=-1, idx7=-1
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%value, has_value )
!!!! WRITE LULOOP, add_variable
      case('units')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%units
!!!!        varname=units, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%units, has_units )
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
end subroutine read_xml_type_density_xml
!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUDEFLT, add_typedef
subroutine init_xml_type_density_xml_array( dvar )
   type(density_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(1) )
end subroutine init_xml_type_density_xml_array
subroutine init_xml_type_density_xml(dvar)
   type(density_xml) :: dvar
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%value = 0.0
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%units = 'atom/b-cm'
end subroutine init_xml_type_density_xml
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITE, add_typedef
subroutine write_xml_type_density_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(density_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_density_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_density_xml_array

subroutine write_xml_type_density_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(density_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_double( info, 'value', indent+3, dvar%value)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'units', indent+3, dvar%units)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_density_xml

!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_nuclide_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar, count_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(nuclide_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar
   integer, intent(inout)                       :: count_dvar

   integer                                      :: newsize
   integer                                      :: size_dvar
   type(nuclide_xml), dimension(:), pointer :: newvar

   count_dvar = count_dvar + 1
   size_dvar = size(dvar)
   if (count_dvar .gt. size_dvar) then
       newsize = size_dvar * 2
       allocate(newvar(1:newsize))
       newvar(1:size_dvar) = dvar
       deallocate(dvar)
       dvar => newvar
   endif

   call read_xml_type_nuclide_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(count_dvar), has_dvar )
end subroutine read_xml_type_nuclide_xml_array

!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_nuclide_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(nuclide_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_name
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_xs
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_ao
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_wo
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable
   has_name                             = .false.
!!!! WRITE LUSTART, add_variable
   has_xs                               = .false.
!!!! WRITE LUSTART, add_variable
   has_ao                               = .false.
!!!! WRITE LUSTART, add_variable
   has_wo                               = .false.
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LULOOP, add_begin_loop
   call init_xml_type_nuclide_xml(dvar)
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
      case('name')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%name
!!!!        varname=name, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%name, has_name )
!!!! WRITE LULOOP, add_variable
      case('xs')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%xs
!!!!        varname=xs, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%xs, has_xs )
!!!! WRITE LULOOP, add_variable
      case('ao')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%ao
!!!!        varname=ao, vartype=double, idx3=13, idx6=-1, idx7=-1
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%ao, has_ao )
!!!! WRITE LULOOP, add_variable
      case('wo')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%wo
!!!!        varname=wo, vartype=double, idx3=13, idx6=-1, idx7=-1
         call read_xml_double( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%wo, has_wo )
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
end subroutine read_xml_type_nuclide_xml
!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUDEFLT, add_typedef
subroutine init_xml_type_nuclide_xml_array( dvar )
   type(nuclide_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(1) )
end subroutine init_xml_type_nuclide_xml_array
subroutine init_xml_type_nuclide_xml(dvar)
   type(nuclide_xml) :: dvar
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%name = ''
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%xs = ''
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%ao = 0.0
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%wo = 0.0
end subroutine init_xml_type_nuclide_xml
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITE, add_typedef
subroutine write_xml_type_nuclide_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(nuclide_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_nuclide_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_nuclide_xml_array

subroutine write_xml_type_nuclide_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(nuclide_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'name', indent+3, dvar%name)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'xs', indent+3, dvar%xs)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_double( info, 'ao', indent+3, dvar%ao)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_double( info, 'wo', indent+3, dvar%wo)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_nuclide_xml

!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_sab_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar, count_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(sab_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar
   integer, intent(inout)                       :: count_dvar

   integer                                      :: newsize
   integer                                      :: size_dvar
   type(sab_xml), dimension(:), pointer :: newvar

   count_dvar = count_dvar + 1
   size_dvar = size(dvar)
   if (count_dvar .gt. size_dvar) then
       newsize = size_dvar * 2
       allocate(newvar(1:newsize))
       newvar(1:size_dvar) = dvar
       deallocate(dvar)
       dvar => newvar
   endif

   call read_xml_type_sab_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(count_dvar), has_dvar )
end subroutine read_xml_type_sab_xml_array

!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_sab_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(sab_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_name
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_xs
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable
   has_name                             = .false.
!!!! WRITE LUSTART, add_variable
   has_xs                               = .false.
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LULOOP, add_begin_loop
   call init_xml_type_sab_xml(dvar)
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
      case('name')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%name
!!!!        varname=name, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%name, has_name )
!!!! WRITE LULOOP, add_variable
      case('xs')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%xs
!!!!        varname=xs, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%xs, has_xs )
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
   if ( .not. has_name ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on name')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_xs ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on xs')
   endif
end subroutine read_xml_type_sab_xml
!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUDEFLT, add_typedef
subroutine init_xml_type_sab_xml_array( dvar )
   type(sab_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(1) )
end subroutine init_xml_type_sab_xml_array
subroutine init_xml_type_sab_xml(dvar)
   type(sab_xml) :: dvar
end subroutine init_xml_type_sab_xml
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITE, add_typedef
subroutine write_xml_type_sab_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(sab_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_sab_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_sab_xml_array

subroutine write_xml_type_sab_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(sab_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'name', indent+3, dvar%name)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'xs', indent+3, dvar%xs)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_sab_xml

!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_material_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar, count_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(material_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar
   integer, intent(inout)                       :: count_dvar

   integer                                      :: newsize
   integer                                      :: size_dvar
   type(material_xml), dimension(:), pointer :: newvar

   count_dvar = count_dvar + 1
   size_dvar = size(dvar)
   if (count_dvar .gt. size_dvar) then
       newsize = size_dvar * 2
       allocate(newvar(1:newsize))
       newvar(1:size_dvar) = dvar
       deallocate(dvar)
       dvar => newvar
   endif

   call read_xml_type_material_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(count_dvar), has_dvar )
end subroutine read_xml_type_material_xml_array

!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_material_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(material_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_id
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_density
!!!! WRITE LUPROLOG, add_variable, idx7>=1
   type(nuclide_xml), dimension(:), pointer :: temp_nuclides => null()
   integer :: count_nuclides
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_nuclides
!!!! WRITE LUPROLOG, add_variable, idx7>=1
   type(nuclide_xml), dimension(:), pointer :: temp_elements => null()
   integer :: count_elements
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_elements
!!!! WRITE LUPROLOG, add_variable, idx7>=1
   type(sab_xml), dimension(:), pointer :: temp_sab => null()
   integer :: count_sab
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_sab
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable
   has_id                               = .false.
!!!! WRITE LUSTART, add_variable
   has_density                          = .false.
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_nuclides = 0
!!!! WRITE LUSTART, add_variable
   has_nuclides                         = .false.
!!!! WRITE LUSTART, add_variable
   allocate(dvar%nuclides(1))
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_elements = 0
!!!! WRITE LUSTART, add_variable
   has_elements                         = .false.
!!!! WRITE LUSTART, add_variable
   allocate(dvar%elements(1))
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_sab = 0
!!!! WRITE LUSTART, add_variable
   has_sab                              = .false.
!!!! WRITE LUSTART, add_variable
   allocate(dvar%sab(1))
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LULOOP, add_begin_loop
   call init_xml_type_material_xml(dvar)
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
      case('density')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%density
!!!!        varname=density, vartype=density_xml, idx3=28, idx6=-1, idx7=-1
         call read_xml_type_density_xml( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%density, has_density )
!!!! WRITE LULOOP, add_variable
      case('nuclide')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%nuclides
!!!!        varname=nuclides, vartype=nuclide_xml-array, idx3=32, idx6=-1, idx7= 4
         call read_xml_type_nuclide_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%nuclides, has_nuclides, count_nuclides )
!!!! WRITE LULOOP, add_variable
      case('element')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%elements
!!!!        varname=elements, vartype=nuclide_xml-array, idx3=32, idx6=-1, idx7= 4
         call read_xml_type_nuclide_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%elements, has_elements, count_elements )
!!!! WRITE LULOOP, add_variable
      case('sab')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%sab
!!!!        varname=sab, vartype=sab_xml-array, idx3=35, idx6=-1, idx7= 3
         call read_xml_type_sab_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%sab, has_sab, count_sab )
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
   if ( .not. has_density ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on density')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_nuclides ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on nuclides')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_nuclides(1:count_nuclides))
   temp_nuclides = dvar%nuclides(1:count_nuclides)
   deallocate(dvar%nuclides)
   dvar%nuclides => temp_nuclides
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_elements ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on elements')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_elements(1:count_elements))
   temp_elements = dvar%elements(1:count_elements)
   deallocate(dvar%elements)
   dvar%elements => temp_elements
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_sab ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on sab')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_sab(1:count_sab))
   temp_sab = dvar%sab(1:count_sab)
   deallocate(dvar%sab)
   dvar%sab => temp_sab
end subroutine read_xml_type_material_xml
!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUDEFLT, add_typedef
subroutine init_xml_type_material_xml_array( dvar )
   type(material_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(1) )
end subroutine init_xml_type_material_xml_array
subroutine init_xml_type_material_xml(dvar)
   type(material_xml) :: dvar
end subroutine init_xml_type_material_xml
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITE, add_typedef
subroutine write_xml_type_material_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(material_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_material_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_material_xml_array

subroutine write_xml_type_material_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(material_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_integer( info, 'id', indent+3, dvar%id)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_xml_type_density_xml( info, 'density', indent+3, dvar%density)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_xml_type_nuclide_xml_array( info, 'nuclide', indent+3, dvar%nuclides)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_xml_type_nuclide_xml_array( info, 'element', indent+3, dvar%elements)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_xml_type_sab_xml_array( info, 'sab', indent+3, dvar%sab)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_material_xml

!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, write_prolog
subroutine read_xml_file_materials_t(fname, lurep, errout)
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
   character(len=1000), dimension(1:1000)  :: data
   integer                                :: nodata
!!!! WRITE LUPROLOG, add_variable, idx7>=1
   type(material_xml), dimension(:), pointer :: temp_material_ => null()
   integer :: count_material_
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_material_
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_default_xs_
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_material_ = 0
!!!! WRITE LUSTART, add_variable
   has_material_                        = .false.
!!!! WRITE LUSTART, add_variable
   allocate(material_(1))
!!!! WRITE LUSTART, add_variable
   has_default_xs_                      = .false.
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LUMAIN, main
!!!! WRITE LUMAIN, write_prolog

   call init_xml_file_materials_t
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
   if ( starttag /= "materials" ) then
      call xml_report_errors( info, &
         'XML-file should have root element "materials"')
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
      case('material')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=material_
!!!!        varname=material_, vartype=material_xml-array, idx3=38, idx6=-1, idx7= 4
         call read_xml_type_material_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            material_, has_material_, count_material_ )
!!!! WRITE LULOOP, add_variable
      case('default_xs')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=default_xs_
!!!!        varname=default_xs_, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            default_xs_, has_default_xs_ )
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
   if ( .not. has_material_ ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=false
      error = .true.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on material_')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_material_(1:count_material_))
   temp_material_ = material_(1:count_material_)
   deallocate(material_)
   material_ => temp_material_
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_default_xs_ ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=false
      error = .true.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on default_xs_')
   endif
   if ( present(errout) ) errout = error
   call xml_close(info)
end subroutine

!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITP, main
!!!! WRITE LUWRV, main
!!!! WRITE LUWRITP, write_prolog
subroutine write_xml_file_materials_t(fname, lurep)
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
      '<materials>'
!!!! WRITE LUWRV, add_variable, component=false
   call write_xml_type_material_xml_array( info, 'material', indent+3, material_)
!!!! WRITE LUWRV, add_variable, component=false
   call write_to_xml_word( info, 'default_xs', indent+3, default_xs_)
   write(info%lun,'(a)') '</materials>'
   call xml_close(info)
end subroutine

!!!! WRITE LUINIT, main
subroutine init_xml_file_materials_t

end subroutine

end module
