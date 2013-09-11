!!!! WRITE LUDEF, main
!!!! WRITE LUDEF, write_prolog
module xml_data_tallies_t
   use READ_XML_PRIMITIVES
   use WRITE_XML_PRIMITIVES
   use XMLPARSE
   implicit none
   save
   integer, private :: lurep_
   logical, private :: strict_
!!!! WRITE LUDEF, add_typedef

type mesh_xml
!!!! WRITE LUDEF, add_variable
   integer                                         :: id
!!!! WRITE LUDEF, add_variable
   character(len=12)                                :: type
!!!! WRITE LUDEF, add_variable
   integer, dimension(:), pointer                  :: dimension => null()
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0)), dimension(:), pointer   :: lower_left => null()
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0)), dimension(:), pointer   :: upper_right => null()
!!!! WRITE LUDEF, add_variable
   real(kind=kind(1.0d0)), dimension(:), pointer   :: width => null()
end type mesh_xml
!!!! WRITE LUDEF, add_typedef

type filter_xml
!!!! WRITE LUDEF, add_variable
   character(len=20)                                :: type
!!!! WRITE LUDEF, add_variable
   character(len=20), dimension(:), pointer         :: bins => null()
!!!! WRITE LUDEF, add_variable
   character(len=20)                                :: groups
end type filter_xml
!!!! WRITE LUDEF, add_typedef

type tally_xml
!!!! WRITE LUDEF, add_variable
   integer                                         :: id
!!!! WRITE LUDEF, add_variable
   character(len=52)                                :: label
!!!! WRITE LUDEF, add_variable
   character(len=12)                                :: estimator
!!!! WRITE LUDEF, add_variable
   type(filter_xml), dimension(:), pointer         :: filter => null()
!!!! WRITE LUDEF, add_variable
   character(len=20), dimension(:), pointer         :: scores => null()
!!!! WRITE LUDEF, add_variable
   character(len=12), dimension(:), pointer         :: nuclides => null()
!!!! WRITE LUDEF, add_variable
   type(filter_xml), dimension(:), pointer         :: filters => null()
end type tally_xml
!!!! WRITE LUDEF, add_variable
   type(mesh_xml), dimension(:), pointer           :: mesh_ => null()
!!!! WRITE LUDEF, add_variable
   type(tally_xml), dimension(:), pointer          :: tally_ => null()
!!!! WRITE LUDEF, add_variable
   character(len=3)                                :: separate_
!!!! WRITE LUSUBS, main
!!!! WRITE LUSUBS, write_prolog
contains
!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_mesh_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar, count_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(mesh_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar
   integer, intent(inout)                       :: count_dvar

   type(mesh_xml), dimension(:), pointer :: temp_dvar

   count_dvar = count_dvar + 1
   do while (count_dvar .gt. size(dvar))
       allocate(temp_dvar(1:size(dvar)*2))
       temp_dvar(1:size(dvar)) = dvar
       deallocate(dvar)
       dvar => temp_dvar
       temp_dvar => null()
   enddo

   call read_xml_type_mesh_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(count_dvar), has_dvar )
end subroutine read_xml_type_mesh_xml_array

!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_mesh_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(mesh_xml), intent(inout)  :: dvar
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
   logical                                         :: has_upper_right
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_width
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
   has_upper_right                      = .false.
!!!! WRITE LUSTART, add_variable
   has_width                            = .false.
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LULOOP, add_begin_loop
   call init_xml_type_mesh_xml(dvar)
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
      case('upper_right')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%upper_right
!!!!        varname=upper_right, vartype=double-array, idx3=15, idx6=-1, idx7=-1
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%upper_right, has_upper_right )
!!!! WRITE LULOOP, add_variable
      case('width')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%width
!!!!        varname=width, vartype=double-array, idx3=15, idx6=-1, idx7=-1
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%width, has_width )
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
   if ( .not. has_upper_right ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on upper_right')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_width ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on width')
   endif
end subroutine read_xml_type_mesh_xml
!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUDEFLT, add_typedef
subroutine init_xml_type_mesh_xml_array( dvar )
   type(mesh_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(1) )
end subroutine init_xml_type_mesh_xml_array
subroutine init_xml_type_mesh_xml(dvar)
   type(mesh_xml) :: dvar
end subroutine init_xml_type_mesh_xml
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITE, add_typedef
subroutine write_xml_type_mesh_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(mesh_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_mesh_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_mesh_xml_array

subroutine write_xml_type_mesh_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(mesh_xml)                      :: dvar
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
   if(associated(dvar%upper_right)) call write_to_xml_double_array( info, 'upper_right', indent+3, dvar%upper_right)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%width)) call write_to_xml_double_array( info, 'width', indent+3, dvar%width)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_mesh_xml

!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_filter_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar, count_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(filter_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar
   integer, intent(inout)                       :: count_dvar

   type(filter_xml), dimension(:), pointer :: temp_dvar

   count_dvar = count_dvar + 1
   do while (count_dvar .gt. size(dvar))
       allocate(temp_dvar(1:size(dvar)*2))
       temp_dvar(1:size(dvar)) = dvar
       deallocate(dvar)
       dvar => temp_dvar
       temp_dvar => null()
   enddo

   call read_xml_type_filter_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(count_dvar), has_dvar )
end subroutine read_xml_type_filter_xml_array

!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_filter_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(filter_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_type
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_bins
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_groups
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable
   has_type                             = .false.
!!!! WRITE LUSTART, add_variable
   has_bins                             = .false.
!!!! WRITE LUSTART, add_variable
   has_groups                           = .false.
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LULOOP, add_begin_loop
   call init_xml_type_filter_xml(dvar)
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
      case('type')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%type
!!!!        varname=type, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type, has_type )
!!!! WRITE LULOOP, add_variable
      case('bins')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%bins
!!!!        varname=bins, vartype=word-array, idx3=19, idx6=-1, idx7=-1
         call read_xml_word_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%bins, has_bins )
!!!! WRITE LULOOP, add_variable
      case('groups')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%groups
!!!!        varname=groups, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%groups, has_groups )
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
   if ( .not. has_bins ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on bins')
   endif
end subroutine read_xml_type_filter_xml
!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUDEFLT, add_typedef
subroutine init_xml_type_filter_xml_array( dvar )
   type(filter_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(1) )
end subroutine init_xml_type_filter_xml_array
subroutine init_xml_type_filter_xml(dvar)
   type(filter_xml) :: dvar
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%type = ''
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%groups = ''
end subroutine init_xml_type_filter_xml
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITE, add_typedef
subroutine write_xml_type_filter_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(filter_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_filter_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_filter_xml_array

subroutine write_xml_type_filter_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(filter_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'type', indent+3, dvar%type)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%bins)) call write_to_xml_word_array( info, 'bins', indent+3, dvar%bins)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'groups', indent+3, dvar%groups)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_filter_xml

!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_tally_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar, count_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(tally_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar
   integer, intent(inout)                       :: count_dvar

   type(tally_xml), dimension(:), pointer :: temp_dvar

   count_dvar = count_dvar + 1
   do while (count_dvar .gt. size(dvar))
       allocate(temp_dvar(1:size(dvar)*2))
       temp_dvar(1:size(dvar)) = dvar
       deallocate(dvar)
       dvar => temp_dvar
       temp_dvar => null()
   enddo

   call read_xml_type_tally_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(count_dvar), has_dvar )
end subroutine read_xml_type_tally_xml_array

!!!! WRITE LUPROLOG, add_typedef
subroutine read_xml_type_tally_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(tally_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_id
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_label
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_estimator
!!!! WRITE LUPROLOG, add_variable, idx7>=1
   type(filter_xml), dimension(:), pointer :: temp_filter => null()
   integer :: count_filter
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_filter
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_scores
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_nuclides
!!!! WRITE LUPROLOG, add_variable, idx7>=1
   type(filter_xml), dimension(:), pointer :: temp_filters => null()
   integer :: count_filters
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_filters
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable
   has_id                               = .false.
!!!! WRITE LUSTART, add_variable
   has_label                            = .false.
!!!! WRITE LUSTART, add_variable
   has_estimator                        = .false.
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_filter = 0
!!!! WRITE LUSTART, add_variable
   has_filter                           = .false.
!!!! WRITE LUSTART, add_variable
   allocate(dvar%filter(1))
!!!! WRITE LUSTART, add_variable
   has_scores                           = .false.
!!!! WRITE LUSTART, add_variable
   has_nuclides                         = .false.
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_filters = 0
!!!! WRITE LUSTART, add_variable
   has_filters                          = .false.
!!!! WRITE LUSTART, add_variable
   allocate(dvar%filters(1))
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LULOOP, add_begin_loop
   call init_xml_type_tally_xml(dvar)
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
      case('label')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%label
!!!!        varname=label, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%label, has_label )
!!!! WRITE LULOOP, add_variable
      case('estimator')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%estimator
!!!!        varname=estimator, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%estimator, has_estimator )
!!!! WRITE LULOOP, add_variable
      case('filter')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%filter
!!!!        varname=filter, vartype=filter_xml-array, idx3=32, idx6=-1, idx7= 3
         call read_xml_type_filter_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%filter, has_filter, count_filter )
!!!! WRITE LULOOP, add_variable
      case('scores')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%scores
!!!!        varname=scores, vartype=word-array, idx3=19, idx6=-1, idx7=-1
         call read_xml_word_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scores, has_scores )
!!!! WRITE LULOOP, add_variable
      case('nuclides')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%nuclides
!!!!        varname=nuclides, vartype=word-array, idx3=19, idx6=-1, idx7=-1
         call read_xml_word_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%nuclides, has_nuclides )
!!!! WRITE LULOOP, add_variable
      case('filters')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=dvar%filters
!!!!        varname=filters, vartype=filter_xml-array, idx3=32, idx6=-1, idx7= 3
         call read_xml_type_filter_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%filters, has_filters, count_filters )
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
   if ( .not. has_filter ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on filter')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_filter(1:count_filter))
   temp_filter = dvar%filter(1:count_filter)
   deallocate(dvar%filter)
   dvar%filter => temp_filter
   temp_filter => null()
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_scores ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on scores')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_nuclides ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on nuclides')
   endif
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_filters ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=true
      has_dvar = .false.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on filters')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_filters(1:count_filters))
   temp_filters = dvar%filters(1:count_filters)
   deallocate(dvar%filters)
   dvar%filters => temp_filters
   temp_filters => null()
end subroutine read_xml_type_tally_xml
!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUDEFLT, add_typedef
subroutine init_xml_type_tally_xml_array( dvar )
   type(tally_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(1) )
end subroutine init_xml_type_tally_xml_array
subroutine init_xml_type_tally_xml(dvar)
   type(tally_xml) :: dvar
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%label = ''
!!!! WRITE LUDFLT, add_variable, idx4>0, compoment=true
   dvar%estimator = ''
end subroutine init_xml_type_tally_xml
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITE, add_typedef
subroutine write_xml_type_tally_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(tally_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_tally_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_tally_xml_array

subroutine write_xml_type_tally_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(tally_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_integer( info, 'id', indent+3, dvar%id)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'label', indent+3, dvar%label)
!!!! WRITE LUWRITE, add_variable, component=true
   call write_to_xml_word( info, 'estimator', indent+3, dvar%estimator)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%filter)) call write_xml_type_filter_xml_array( info, 'filter', indent+3, dvar%filter)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%scores)) call write_to_xml_word_array( info, 'scores', indent+3, dvar%scores)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%nuclides)) call write_to_xml_word_array( info, 'nuclides', indent+3, dvar%nuclides)
!!!! WRITE LUWRITE, add_variable, component=true
   if(associated(dvar%filters)) call write_xml_type_filter_xml_array( info, 'filters', indent+3, dvar%filters)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_tally_xml

!!!! WRITE LUPROLOG, open_tmp_files
!!!! WRITE LUPROLOG, write_prolog
subroutine read_xml_file_tallies_t(fname, lurep, errout)
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
   type(mesh_xml), dimension(:), pointer :: temp_mesh_ => null()
   integer :: count_mesh_
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_mesh_
!!!! WRITE LUPROLOG, add_variable, idx7>=1
   type(tally_xml), dimension(:), pointer :: temp_tally_ => null()
   integer :: count_tally_
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_tally_
!!!! WRITE LUPROLOG, add_variable
   logical                                         :: has_separate_
!!!! WRITE LUSTART, open_tmp_files
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_mesh_ = 0
!!!! WRITE LUSTART, add_variable
   has_mesh_                            = .false.
!!!! WRITE LUSTART, add_variable
   allocate(mesh_(1))
!!!! WRITE LUSTART, add_variable, idx7>=1
   count_tally_ = 0
!!!! WRITE LUSTART, add_variable
   has_tally_                           = .false.
!!!! WRITE LUSTART, add_variable
   allocate(tally_(1))
!!!! WRITE LUSTART, add_variable
   has_separate_                        = .false.
!!!! WRITE LULOOP, open_tmp_files
!!!! WRITE LUMAIN, main
!!!! WRITE LUMAIN, write_prolog

   call init_xml_file_tallies_t
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
   if ( starttag /= "tallies" ) then
      call xml_report_errors( info, &
         'XML-file should have root element "tallies"')
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
      case('mesh')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=mesh_
!!!!        varname=mesh_, vartype=mesh_xml-array, idx3=29, idx6=-1, idx7= 4
         call read_xml_type_mesh_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            mesh_, has_mesh_, count_mesh_ )
!!!! WRITE LULOOP, add_variable
      case('tally')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=tally_
!!!!        varname=tally_, vartype=tally_xml-array, idx3=35, idx6=-1, idx7= 4
         call read_xml_type_tally_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            tally_, has_tally_, count_tally_ )
!!!! WRITE LULOOP, add_variable
      case('assume_separate')
!!!! WRITE LULOOP, add_variable, idx6<=0, varcomp=separate_
!!!!        varname=separate_, vartype=word, idx3=17, idx6=-1, idx7=-1
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            separate_, has_separate_ )
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
   if ( .not. has_mesh_ ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=false
      error = .true.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on mesh_')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_mesh_(1:count_mesh_))
   temp_mesh_ = mesh_(1:count_mesh_)
   deallocate(mesh_)
   mesh_ => temp_mesh_
   temp_mesh_ => null()
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_tally_ ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=false
      error = .true.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on tally_')
   endif
!!!! WRITE LUEND, add_variable, idx4<=0, idx7>=0
   allocate(temp_tally_(1:count_tally_))
   temp_tally_ = tally_(1:count_tally_)
   deallocate(tally_)
   tally_ => temp_tally_
   temp_tally_ => null()
!!!! WRITE LUEND, add_variable, idx4 <= 0
   if ( .not. has_separate_ ) then
!!!! WRITE LUEND, add_variable, idx4<=0, component=false
      error = .true.
!!!! WRITE LUEND, add_variable, idx4<=0
      call xml_report_errors(info, 'Missing data on separate_')
   endif
   if ( present(errout) ) errout = error
   call xml_close(info)
end subroutine

!!!! WRITE LUDEFLT, open_tmp_files
!!!! WRITE LUWRITE, open_tmp_files
!!!! WRITE LUWRITP, main
!!!! WRITE LUWRV, main
!!!! WRITE LUWRITP, write_prolog
subroutine write_xml_file_tallies_t(fname, lurep)
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
      '<tallies>'
!!!! WRITE LUWRV, add_variable, component=false
   if(associated(mesh_))   call write_xml_type_mesh_xml_array( info, 'mesh', indent+3, mesh_)
!!!! WRITE LUWRV, add_variable, component=false
   if(associated(tally_))   call write_xml_type_tally_xml_array( info, 'tally', indent+3, tally_)
!!!! WRITE LUWRV, add_variable, component=false
   call write_to_xml_word( info, 'assume_separate', indent+3, separate_)
   write(info%lun,'(a)') '</tallies>'
   call xml_close(info)
end subroutine

!!!! WRITE LUINIT, main
subroutine init_xml_file_tallies_t

end subroutine

end module
